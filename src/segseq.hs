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

processSites :: Handle -> Handle -> ([String], [String])-> IO()
processSites fh rh (f,r)= do
                             forkIO $ do mapM_ (hPutStrLn fh) f
                                         hClose fh
                             forkIO $ do mapM_ (hPutStrLn rh) r
                                         hClose rh
                             return ()


run :: Settings -> [Annotation] -> IO ()
run s a=   do
             createIndexFromFasta (fasta s)
             (sinF, soutF, serrF, pidF) <- runInteractiveCommand   $ " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d " ++ (fasta s) ++ " -R "
             (sinR, soutR, serrR, pidR) <- runInteractiveCommand   $ " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d " ++ (fasta s) ++ " -R | reverse_complement.pl "
             hSetBinaryMode sinF False
             hSetBinaryMode sinR False
             return  (extractContent s a) >>= (processSites sinF sinR)
             hGetContents soutF >>= putStr
             hGetContents soutR >>= putStr
             hGetContents serrF >>= hPutStrLn stderr
             hGetContents serrR >>= hPutStrLn stderr
             waitForProcess pidF
             waitForProcess pidR
             return ();
