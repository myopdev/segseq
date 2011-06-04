module FastaDB (indexFastaFile, getSequence) where

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString.Lazy.Char8 as C
import Database.HDBC
import Database.HDBC.Sqlite3
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
import System.IO
import System.IO.Unsafe
import Data.List (isInfixOf, isPrefixOf, foldl')
import Data.List.Split
import Foreign


trim      :: String -> String
trim      = f . f
   where f = reverse . dropWhile isSpace

getSequence :: FilePath -> FilePath -> String -> IO(L.ByteString)
getSequence db fasta code =
  do words <- return(splitOn " " (trim $! cleanSpaces  code))
     seqname <- return(head words)
     positions <- return $! (getPositions (tail words))
     if (length(positions) <= 0)
        then do return (L.empty)
        else do seqByPositions db fasta seqname positions
  where getPositions [] = []
        getPositions [x] = []
        getPositions [x,y] = [((read x)::Int64, (read y ):: Int64)]
        getPositions (x:y:rest) = [((read x)::Int64, (read y) :: Int64)] ++ (getPositions rest)


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
     db <- connectSqlite3 absPathDB
     begin <- getInt db key $ getStartPositionFromDatabase
     small <- getInt db key $ getSmallestOffsetFromDatabase
     largest <- getInt db key $ getLargestOffsetFromDatabase
     beginOfSequence <- return(begin + small + 1)
     hdl <- openFile absPathFasta ReadMode
     hSeek hdl AbsoluteSeek (fromIntegral beginOfSequence)
     fsize <- hFileSize hdl
     seqdata <- if ((largest - small - 1) > 0 )
                    then do fp <- mallocForeignPtrBytes (fromIntegral (largest - small -1))
                            len <- withForeignPtr fp $ \buf -> hGetBuf hdl buf (fromIntegral (largest - small - 1))
                            lazySlurp fp (fromIntegral 0)  (fromIntegral len)
                    else do return( L.empty)
     hClose hdl
     disconnect db
     return (seqdata)

sumOffset :: Connection -> String -> Int64 -> Int64 -> IO (Int64)
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
     db <- connectSqlite3 absPathDB
     begin <- getInt db key $ getStartPositionFromDatabase
     small <- getInt db key $ getSmallestOffsetFromDatabase
     largest <- getInt db key $ getLargestOffsetFromDatabase
     beginOfSequence <- return(begin + small + 1)
     endOfSequence <- return(begin + largest )
     out <- if (beginOfSequence <= endOfSequence) 
               then do offsetStart <- sumOffset db key b 0
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
                       seqdata <- if ((subSeqEnd - subSeqBegin + 1) > 0)
                                     then do fp <- mallocForeignPtrBytes (fromIntegral (subSeqEnd - subSeqBegin + 1))
                                             len <- withForeignPtr fp $ \buf -> hGetBuf hdl buf (fromIntegral (subSeqEnd - subSeqBegin + 1))
                                             lazySlurp fp 0 (fromIntegral len)
                                     else do return(L.empty)
                       hClose hdl
                       return seqdata
               else do hPutStrLn stderr $ "ERROR: could not find sequence: " ++ key ++ ":" ++ (show b) ++ "," ++ (show e)
                       return L.empty 
     disconnect db
     return out
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
       case (isSpace (w2c w)) of
          True ->   loop (len-1) p (acc)
          False ->  loop (len-1) p (w `L.cons` acc)




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
  do a <- connectSqlite3  absPath
     run a "create table if not exists StartPosition ( key string, value int64, primary key (key, value) )" []
     run a "create table if not exists Offset (key string, value int64, type bool,  primary key (key, value) )" []
     commit a
     withTransaction a insertIntoDatabase
     run a "delete from StartPosition where key = \"\"" []
     run a "delete from Offset where key = \"\"" []
     commit a
     disconnect a
     return()
 where insertIntoDatabase a = ((L.readFile absPathFasta) >>= findKeysFromContent a)

findKeysFromContent :: Connection -> L.ByteString -> IO()
findKeysFromContent  a content =
  do endLine <- return(L.pack([c2w '\n']))
     content <-return( L.append content endLine)
     startOfSeq <- return  (L.findIndices (\x -> ((c2w '>')  == x)) content)
     seqEntries <- return  (tail (L.split (c2w '>') content) )
     newLines <- return (map ( L.findIndices (\x ->  ((x == c2w '\r') || (x == c2w '\n'))   )) seqEntries)
     spaces <- return (map ( L.findIndices (\x ->  ((x == c2w '\t') || (x == c2w ' '))   )) seqEntries)
     seqNames <- return  (map (L.takeWhile (\x -> (not ((x == c2w ' ') || (x == c2w '\r') || (x == c2w '\n')))) ) seqEntries)
     startPositionByName <- return $! (zip seqNames startOfSeq)
     newLinesByName <- return  $! (zip seqNames newLines)
     spacesByName <- return $! (zip seqNames spaces)
     forM startPositionByName (addKeyStartPosition a)
     forM newLinesByName (addKeysNewLineList a)
     forM spacesByName (addKeysSpaceList a)
     return ()

addKeyStartPosition :: Connection -> (L.ByteString, Int64)  -> IO()
addKeyStartPosition a (name, pos) =
  do run a "INSERT INTO StartPosition (key, value) VALUES (?,?)" [toSql $ C.unpack (name), toSql pos]
     return ()

addKeyNewLine :: Connection -> (L.ByteString, Int64) -> IO()
addKeyNewLine a (name, pos) =
  do x <- return (pos + 1)
     handleSqlError $ run a "INSERT INTO Offset (key, value, type) VALUES (?,?,?) " [toSql $ C.unpack (name), (toSql x), (toSql True)]
     return ()

addKeySpace :: Connection -> (L.ByteString, Int64) -> IO()
addKeySpace a (name, pos) =
  do run a "INSERT INTO Offset (key, value, type) VALUES (?,?,?)" [toSql $ C.unpack (name), toSql  (pos + 1), toSql False]
     return ()


addKeysNewLineList ::Connection -> (L.ByteString, [Int64]) -> IO()
addKeysNewLineList a (name, list) =
  do forM list (\ x -> addKeyNewLine a (name, x))
     return()

addKeysSpaceList ::Connection -> (L.ByteString, [Int64]) -> IO()
addKeysSpaceList a (name, list) =
  do forM list (\ x -> addKeySpace a (name, x))
     return()



getKeys :: Connection -> IO [String]
getKeys db =
  do r <- getKeysFromDatabase db
     return (map (\tuple -> fromSql (tuple!!0)::String ) r)


getInt :: Connection -> String -> (Connection -> String -> IO ([[SqlValue]]) ) ->IO (Int64)
getInt db key k =
     do x <- k db key
        if ((length x) > 0) 
           then case (x!!0!!0) of 
                     SqlNull -> return (0);
                     _ -> return ((fromSql (x!!0!!0))::Int64)
           else return 0


getKeysFromDatabase  :: Connection -> IO ([[SqlValue]])
getKeysFromDatabase db =
  do quickQuery' db ("SELECT DISTINCT key FROM StartPosition") []

getStartPositionFromDatabase  :: Connection -> String -> IO ([[SqlValue]])
getStartPositionFromDatabase db key =
  do quickQuery' db ("SELECT value FROM StartPosition where key = ?") [toSql key]

getLargestOffsetFromDatabase  :: Connection -> String -> IO ([[SqlValue]])
getLargestOffsetFromDatabase db key =
  do quickQuery' db ("SELECT max(value) FROM Offset where key =? and type= 'True'") [toSql key]

getSmallestOffsetFromDatabase  :: Connection -> String -> IO ([[SqlValue]])
getSmallestOffsetFromDatabase db key =
  do quickQuery' db ("SELECT min(value) FROM Offset where key = ? and type = 'True'") [toSql key]

getOffsetCountFromDatabase  :: Int64 -> Connection -> String  -> IO ([[SqlValue]])
getOffsetCountFromDatabase pos db key  =
  do quickQuery' db ("SELECT count(value) FROM Offset  where key = ?   and value <= ? and value > (SELECT (min(value)) FROM Offset where key = ? and type = 'True')") [toSql key, SqlInt64 pos, toSql key]


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


