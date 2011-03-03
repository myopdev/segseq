module Main where

import System.IO
import System.Process
import System.Environment
import SequenceAnnotation
import System.Console.GetOpt
import SegSequence
import System.Cmd
import Control.Monad
import GTF
import Data.Maybe
import System.Directory(getTemporaryDirectory, removeFile)


main :: IO()
main =  do
           args <- getArgs
           s <- (liftM (buildSettings (Settings "" "" "" 9 3))) $ (compilerOpts args)
           a <- (liftM (createAnnotationList . readCSV)) $ readFile (gtf s)
           run s a
           return ()


run :: Settings -> [Annotation] -> IO ()
run s a=   do
             createIndexFromFasta (fasta s)

             (sinF, soutF, serrF, pidF) <- runInteractiveCommand   $ " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d " ++ (fasta s) ++ " -R "
             (sinR, soutR, serrR, pidR) <- runInteractiveCommand   $ " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d " ++ (fasta s) ++ " -R "
             (rin, rout, rerr, rpid) <- runInteractiveCommand $ " reverse_complement.pl "
             extractSite (siteName s) sinF sinR s a
             f <- hGetContents soutF
             r <- hGetContents soutR
             hPutStr rin r
             reversed <- hGetContents rout
             putStr f
             putStr reversed
             waitForProcess pidF
             waitForProcess pidR
             waitForProcess rpid
             return ();

