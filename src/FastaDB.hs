module FastaDB (indexFastaFile, getSequence) where

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString.Lazy.Char8 as C
import Database.SQLite
import Data.ByteString.Internal
import Control.Monad
import Data.Char
import Data.Int
import System.FilePath.Posix
import System.Directory
import System.Posix.Directory
import System.Posix.User
import System.Posix.Process
import System.Directory
-- import System.Time
-- import System.Path
-- import System.Cmd.Utils
import System.IO
import System.IO.Unsafe
import Data.List (isInfixOf, isPrefixOf)
import Data.List.Split
import Foreign


trim      :: String -> String
trim      = f . f
   where f = reverse . dropWhile isSpace

getSequence :: FilePath -> FilePath -> String -> IO(L.ByteString)
getSequence db fasta code =
  do words <- return(splitOn " " (trim $ cleanSpaces  code))
     seqname <- return(head words)
     positions <- return (getPositions (tail words))
     if (length(positions) <= 0)
        then do return (L.empty)
        else do seqByPositions db fasta seqname positions
  where getPositions [] = []
        getPositions [x] = []
        getPositions [x,y] = [((read x)::Int64, (read y ):: Int64)]
        getPositions (x:y:rest) = [((read x)::Int64, (read y) :: Int64)] ++ getPositions rest


cleanSpaces :: String -> String
cleanSpaces [] = []
cleanSpaces [x] | (x == ' ') = []
                | otherwise = [x]
cleanSpaces [x,y] = [x] ++ cleanSpaces [y]
cleanSpaces (x:y:xs) | (x == ' ') && (y == ' ') =  (cleanSpaces (y:xs))
                     | otherwise = [x] ++ (cleanSpaces (y:xs))

seqByPositions :: FilePath -> FilePath -> String -> [(Int64, Int64)] -> IO(L.ByteString)
seqByPositions db fasta seqname [] =  return(L.empty)
seqByPositions db fasta seqname [(p1,p2)] = do getSubSequence' db fasta seqname p1 p2
seqByPositions db fasta seqname ((p1,p2):rest) =
  do a<- getSubSequence' db fasta seqname p1 p2
     b<- seqByPositions db fasta seqname rest
     c <- return(L.append a b)
     return (c)

getSequence' :: FilePath -> FilePath -> String -> IO(L.ByteString)
getSequence' fasta dbfile key =
  do absPathFasta <- absDir(fasta)
     absPathDB <- absDir (dbfile)
     db <- openConnection absPathDB
     begin <- getInt db key $ getStartPositionFromDatabase
     small <- getInt db key $ getSmallestOffsetFromDatabase
     largest <- getInt db key $ getLargestOffsetFromDatabase
     beginOfSequence <- return(begin + small + 1)
     hdl <- openFile absPathFasta ReadMode
     hSeek hdl AbsoluteSeek (fromIntegral beginOfSequence)
     fsize <- hFileSize hdl
     fp <- mallocForeignPtrBytes (fromIntegral (largest - small -1))
     len <- withForeignPtr fp $ \buf -> hGetBuf hdl buf (fromIntegral (largest - small - 1))
     seqdata <- lazySlurp fp (fromIntegral 0)  (fromIntegral len)
     hClose hdl
     closeConnection db
     return (seqdata)

sumOffset :: SQLiteHandle -> String -> Int64 -> Int64 -> IO (Int64)
sumOffset db key pos o =
  do small <- getInt db key $ getSmallestOffsetFromDatabase
     x <- getInt db key $ getOffsetCountFromDatabase $fromIntegral (pos + small + o)
     if (not (x == o))
        then do sumOffset db key pos (x)
        else do return x



getSubSequence' :: FilePath -> FilePath -> String -> Int64 -> Int64 -> IO(L.ByteString)
getSubSequence' fasta dbfile key b e =
  do absPathFasta <- absDir(fasta)
     absPathDB <- absDir (dbfile)
     db <- openConnection absPathDB
     begin <- getInt db key $ getStartPositionFromDatabase
     small <- getInt db key $ getSmallestOffsetFromDatabase
     largest <- getInt db key $ getLargestOffsetFromDatabase

     beginOfSequence <- return(begin + small + 1)
     endOfSequence <- return(begin + largest )

     offsetStart <- sumOffset db key b 0
     offsetEnd <- sumOffset db key e 0
     subSeqBegin <- return(offsetStart + b + beginOfSequence - 1)
     subSeqEnd <- return(offsetEnd + e + beginOfSequence - 1)
     subSeqBegin <- fixLimitStart subSeqBegin beginOfSequence
     subSeqBegin <- fixLimitEnd subSeqBegin endOfSequence
     subSeqEnd <- fixLimitStart subSeqEnd beginOfSequence
     subSeqEnd <- fixLimitEnd subSeqEnd endOfSequence

     hdl <- openFile absPathFasta ReadMode
     hSeek hdl AbsoluteSeek (fromIntegral subSeqBegin)
     fsize <- hFileSize hdl
     fp <- mallocForeignPtrBytes (fromIntegral (subSeqEnd - subSeqBegin + 1))
     len <- withForeignPtr fp $ \buf -> hGetBuf hdl buf (fromIntegral (subSeqEnd - subSeqBegin + 1))
     seqdata <- lazySlurp fp 0 (fromIntegral len)
     hClose hdl
     closeConnection db
     return seqdata
  where fixLimitStart pos1 start  =
         do if (pos1 < start)
               then return (start)
               else return pos1
        fixLimitEnd pos1 end =
         do if (pos1 > end)
               then return (end)
               else return pos1




buf_size = 4096 :: Int

lazySlurp :: ForeignPtr Word8 -> Int -> Int -> IO L.ByteString
lazySlurp fp ix len
  | fp `seq` False = undefined
  | ix >= len = return L.empty
  | otherwise = do
      cs <- unsafeInterleaveIO (lazySlurp fp (ix + buf_size) len)
      ws <- withForeignPtr fp $ \p -> loop (min (len-ix) buf_size - 1) 
      	    		      ((p :: Ptr Word8) `plusPtr` ix) cs
      return ws
 where
  loop :: Int -> Ptr Word8 -> L.ByteString -> IO L.ByteString
  loop len p acc
    | len `seq` p `seq` False = undefined
    | len < 0 = return acc
    | otherwise = do
       w <- peekElemOff p len
       loop (len-1) p (w `L.cons` acc)




indexFastaFile :: FilePath -> FilePath-> IO()
indexFastaFile path outfile =
  do absPathFasta <- absDir(path)
     absPathDB <- absDir (outfile)
     x <- getModificationTime absPathFasta
     exist <- doesFileExist absPathDB
     if (not (exist))
        then do
         updateIndex  absPathDB absPathFasta
         return()
        else do
         y <- getModificationTime absPathDB
         if ( x > y )
            then do
             removeFile absPathDB
             updateIndex  absPathDB absPathFasta
             return()
            else do return()

updateIndex ::  FilePath  -> FilePath -> IO()
updateIndex  absPath absPathFasta =
  do a <- openConnection  absPath
     execStatement_ a "BEGIN;"
     execStatement_ a "create table if not exists StartPosition ( key string, value int64, primary key (key, value) )"
     execStatement_ a "create table if not exists Offset (key string, value int64, type bool,  primary key (key, value) )"
     (L.readFile absPathFasta) >>= findKeysFromContent a
     execStatement_ a "delete from StartPosition where key = \"\""
     execStatement_ a "delete from Offset where key = \"\""
     execStatement_ a "END"
     closeConnection a
     return()


findKeysFromContent :: SQLiteHandle -> L.ByteString -> IO()
findKeysFromContent  a content =
  do endLine <- return(L.pack([c2w '\n']))
     content <-return( L.append content endLine)
     startOfSeq <- return (L.findIndices (\x -> ((c2w '>')  == x)) content)
     seqEntries <- return(tail (L.split (c2w '>') content) )
     newLines <- return (map ( L.findIndices (\x ->  ((x == c2w '\r') || (x == c2w '\n'))   )) seqEntries)
     spaces <- return (map ( L.findIndices (\x ->  ((x == c2w '\t') || (x == c2w ' '))   )) seqEntries)
     seqNames <- return (map (L.takeWhile (\x -> (not ((x == c2w '\r') || (x == c2w '\n')))) ) seqEntries)
     startPositionByName <- return (zip seqNames startOfSeq)
     newLinesByName <- return (zip seqNames newLines)
     spacesByName <- return (zip seqNames spaces)
     forM startPositionByName (addKeyStartPosition a)
     forM newLinesByName (addKeysNewLineList a)
     forM spacesByName (addKeysSpaceList a)
     return ()

addKeyStartPosition :: SQLiteHandle -> (L.ByteString, Int64)  -> IO()
addKeyStartPosition a (name, pos) =
  do execStatement_ a $ "INSERT INTO StartPosition (key, value) VALUES (\"" ++ C.unpack(name) ++ "\"," ++  (show pos) ++ ")"
     return ()

addKeyNewLine :: SQLiteHandle -> (L.ByteString, Int64) -> IO()
addKeyNewLine a (name, pos) =
  do execStatement_ a $ "INSERT INTO Offset (key, value, type) VALUES (\"" ++ C.unpack(name) ++ "\"," ++  (show (pos+1)) ++ ", 1)"
     return ()

addKeySpace :: SQLiteHandle -> (L.ByteString, Int64) -> IO()
addKeySpace a (name, pos) =
  do execStatement_ a $ "INSERT INTO Offset (key, value, type) VALUES (\"" ++ C.unpack(name) ++ "\"," ++  (show (pos+1)) ++ ", 0)"
     return ()


addKeysNewLineList ::SQLiteHandle -> (L.ByteString, [Int64]) -> IO()
addKeysNewLineList a (name, list) =
  do forM list (\ x -> addKeyNewLine a (name, x))
     return()

addKeysSpaceList ::SQLiteHandle -> (L.ByteString, [Int64]) -> IO()
addKeysSpaceList a (name, list) =
  do forM list (\ x -> addKeySpace a (name, x))
     return()



getKeys :: SQLiteHandle -> IO [String]
getKeys db =
  do r <- getKeysFromDatabase db
     case r of
      Left a ->  return([])
      Right a -> return (map (\tuple -> snd(tuple !! 0) ) (a !! 0))


getInt :: SQLiteHandle -> String -> (SQLiteHandle -> String -> IO (Either String [[Row Value]]) ) ->IO (Int64)
getInt db key k =
  do x <- k db key
     case x of
       Left a -> return (-1)
       Right a -> case (snd(a!!0!!0!!0)) of
                   Int b -> return(b)
                   Null -> return (-1)




getKeysFromDatabase  :: SQLiteHandle -> IO (Either String [[Row String]])
getKeysFromDatabase db =
  do execStatement db "SELECT DISTINCT key FROM StartPosition"

getStartPositionFromDatabase  :: SQLiteHandle -> String -> IO (Either String [[Row Value]])
getStartPositionFromDatabase db key =
  do execStatement db $ "SELECT value FROM StartPosition where key =\"" ++ key ++ "\""

getLargestOffsetFromDatabase  :: SQLiteHandle -> String -> IO (Either String [[Row Value]])
getLargestOffsetFromDatabase db key =
  do execStatement db $ "SELECT max(value) FROM Offset where key =\"" ++ key ++ "\" and type=1"

getSmallestOffsetFromDatabase  :: SQLiteHandle -> String -> IO (Either String [[Row Value]])
getSmallestOffsetFromDatabase db key =
  do execStatement db $ "SELECT min(value) FROM Offset where key =\"" ++ key ++ "\" and type = 1"

getOffsetCountFromDatabase  :: Int64 -> SQLiteHandle -> String  -> IO (Either String [[Row Value]])
getOffsetCountFromDatabase pos db key  =
  do execStatement db $ "SELECT count(value) FROM Offset  where key =\"" ++ key ++ "\" and value <= " ++ (show (pos)) ++ " and value > (SELECT (min(value)) FROM Offset where key = \"" ++ key ++ "\" and type = 1)"


absDir :: String -> IO String
absDir d
        | "~" `isPrefixOf` d = expandt d
        | otherwise = do
                return $  d
        where
                homedir u = (homeDirectory u) ++ "/"
                myhomedir = do
                        uid <- getEffectiveUserID
                        u <- getUserEntryForID uid
                        return $ homedir u
                expandt [] = return ""
                expandt ('~':'/':cs) = do
                        h <- myhomedir
                        return $ h ++ cs
                expandt (c:cs) = do
                        v <- expandt cs
                        return (c:v)
                findname n [] = (n, "")
                findname n (c:cs)
                         | c == '/' = (n, cs)
                         | otherwise = findname (n++[c]) cs


