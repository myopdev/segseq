module Main where

import System.IO
import System.Process
import System.Environment
import SequenceAnnotation
import System.Console.GetOpt
import SegSequence
import System.Cmd
import Control.Monad
import Control.DeepSeq
import GTF
import Data.Maybe
import System.Directory(getTemporaryDirectory, removeFile)
import Control.Concurrent
import qualified Data.ByteString.Char8 as B

main :: IO()
main =  do
           args <- getArgs
           s <- (liftM (buildSettings (Settings "" "" "" 9 3))) $ (compilerOpts args)
           gtfh <- openFile (gtf s) ReadMode
           (B.hGetContents gtfh) >>= (return . createAnnotationList . readCSV . B.unpack) >>= (run s)
           hClose gtfh
           return ()


getSequenceFromCDBYank :: Settings -> String -> (B.ByteString -> IO(B.ByteString)) -> String-> IO(B.ByteString)
getSequenceFromCDBYank s strand k key  = do
                              (sin, sout, serr, pid) <- runInteractiveCommand   $ strandSpecificCommand
                              hSetBinaryMode sin False
                              forkIO $ do return(key) >>= hPutStrLn sin
                                          hFlush sin
                                          hClose sin
                              ret <- (B.hGetContents sout)
                              waitForProcess pid
                              k (ret)
                    where strandSpecificCommand = case (strand == "+") of
                                                    True -> " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d  " ++ (fasta s) ++ " -R -x 2> err "
                                                    False -> " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d  " ++ (fasta s) ++ " -R -x 2>err| reverse_complement.pl "



run :: Settings -> [Annotation] -> IO ()
run s a=   do
             createIndexFromFasta (fasta s)
             case ( feature s == "cds" ) of
                  True -> getJoinableFeature s a
                  False -> getSimpleFeature s a
             return ();
getJoinableFeature :: Settings -> [Annotation] -> IO()
getJoinableFeature s a = do
                            listCDS <- return (extractCDS s a)
                            mergeListCDS s listCDS
                            return()

mergeListCDS :: Settings-> [([String], [String])] -> IO()
mergeListCDS s [] = return()
mergeListCDS s ((f,r):rest) = do  a<-getForwardCDS s f
                                  b<-getReverseCDS s r
                                  return ((B.pack . (mergeCDS 0) . lines . B.unpack) a) >>= B.putStr
                                  return ((B.pack . (mergeCDS 0) . lines . B.unpack) b) >>= B.putStr
                                  mergeListCDS s rest
                                  return()
mergeCDS ::Integer -> [String]-> String
mergeCDS c []  = ""
mergeCDS c (s:rest)  = case ( (s !! 0) == '>')   of
                       True -> case (c == 0) of
                                 True -> s ++ "\n" ++ mergeCDS (c+1) rest
                                 False ->  mergeCDS (c+1) rest
                       False -> s ++ "\n" ++ mergeCDS c rest

getForwardCDS :: Settings -> [String] -> IO(B.ByteString)
getForwardCDS s [] = return(B.pack "")
getForwardCDS s (l:rest) =  do
                           a <- (getSequenceFromCDBYank s "+" (\x -> return (x)))  l
                           b <- getForwardCDS s rest
                           return (a `B.append` b)

getReverseCDS :: Settings -> [String] -> IO(B.ByteString)
getReverseCDS s [] = return(B.pack "")
getReverseCDS s (l:rest) =  do
                           a <- (getSequenceFromCDBYank s "-" (\x -> return (x))) l
                           b <- getReverseCDS s rest
                           return (b `B.append` a)

getSimpleFeature :: Settings -> [Annotation] -> IO()
getSimpleFeature s a = do
                 (f,r) <- return(extractContent s a)
                 printSimpleFeature s "+" f
                 printSimpleFeature s "-" r
                 return()

printSimpleFeature :: Settings -> String -> [String] -> IO()
printSimpleFeature s strand [] = return()
printSimpleFeature s strand (x:xs) = do
                                        (getSequenceFromCDBYank s strand printAndReturn) x
                                        printSimpleFeature s strand xs
                                        return()

printAndReturn :: B.ByteString -> IO (B.ByteString)
printAndReturn s = do
                      B.putStr s
                      return(s)
