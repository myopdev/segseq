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
import Directory
main :: IO()
main =  do
           s <-  (getArgs >>= (return . compilerOpts) >>= (liftM (buildSettings (Settings "" "" "" 9 3 (-1)))))
           gtfh <- openFile (gtf s) ReadMode
           (B.hGetContents gtfh) >>= (return . createAnnotationList . readCSV . B.unpack) >>= (run s)
           hClose gtfh
           return ()


printFeature :: Settings -> ([String],[String]) -> IO()
printFeature s (f,r)  = do
                           (sinF, soutF, serrF, pidF) <- runInteractiveCommand $ "segseq-getseq -f " ++ (fasta s)
                           (sinR, soutR, serrR, pidR) <- runInteractiveCommand $ "segseq-getseq -f " ++ (fasta s) ++ " | segseq-complement "
                           forkIO $ do return (concat f) >>= hPutStr sinF
                                       hClose sinF
                           forkIO $ do return (concat r) >>= hPutStr sinR
                                       hClose sinR
                           hGetContents soutF >>= putStr
                           hGetContents soutR >>= putStr
                           waitForProcess pidF
                           waitForProcess pidR
                           return()



run :: Settings -> [Annotation] -> IO ()
run s a=   do
             features <- return(extractContent s a)
             printFeature s features
             return ();
