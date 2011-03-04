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

main :: IO()
main =  do
           args <- getArgs
           s <- (liftM (buildSettings (Settings "" "" "" 9 3))) $ (compilerOpts args)
           gtfh <- openFile (gtf s) ReadMode
           (hGetContents gtfh) >>= (return . createAnnotationList . readCSV) >>= (run s)
           hClose gtfh
           return ()


getSequenceFromCDBYank :: Settings -> String -> (String -> IO(String)) -> String-> IO(String)
getSequenceFromCDBYank s strand k key  = do
                              (sin, sout, serr, pid) <- runInteractiveCommand   $ strandSpecificCommand
                              hSetBinaryMode sin False
                              forkIO $ do return(key) >>= hPutStrLn sin
                                          hFlush sin
                                          hClose sin
                              ret <- (hGetContents sout)
                              ret `deepseq` waitForProcess pid
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
                                  return (((mergeCDS 0) . lines) a) >>= putStr
                                  return (((mergeCDS 0) . lines) b) >>= putStr
                                  mergeListCDS s rest
                                  return()
mergeCDS ::Integer -> [String]-> String
mergeCDS c []  = ""
mergeCDS c (s:rest)  = case ( (s !! 0) == '>')   of
                       True -> case (c == 0) of
                                 True -> s ++ "\n" ++ mergeCDS (c+1) rest
                                 False ->  mergeCDS (c+1) rest
                       False -> s ++ "\n" ++ mergeCDS c rest

getForwardCDS :: Settings -> [String] -> IO(String)
getForwardCDS s [] = return("")
getForwardCDS s (l:rest) =  do
                           a <- (getSequenceFromCDBYank s "+" (\x -> return (x)))  l
                           b <- getForwardCDS s rest
                           return (a ++ b)

getReverseCDS :: Settings -> [String] -> IO(String)
getReverseCDS s [] = return("")
getReverseCDS s (l:rest) =  do
                           a <- (getSequenceFromCDBYank s "-" (\x -> return (x))) l
                           b <- getReverseCDS s rest
                           return (b ++ a)

getSimpleFeature :: Settings -> [Annotation] -> IO()
getSimpleFeature s a = do
                 (f,r) <- return(extractContent s a)
                 mapM_ (getSequenceFromCDBYank s "+" printAndReturn) f
                 mapM_ (getSequenceFromCDBYank s "-" printAndReturn) r
                 return()

printAndReturn :: String -> IO (String)
printAndReturn s = do
                      putStr s
                      return(s)
