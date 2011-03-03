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
           s <- (liftM (buildSettings (Settings "" ""  9 3))) $ (compilerOpts args)
           annotationList <- (liftM (createAnnotationList . readCSV)) $ readFile (gtf s)
           createIndexFromFasta (fasta s)

           (sinF, soutF, serrF, pidF) <- runInteractiveCommand   $ " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d " ++ (fasta s) ++ " -R "
           (sinR, soutR, serrR, pidR) <- runInteractiveCommand   $ " cdbyank " ++ (fasta s) ++ ".cidx " ++ " -d " ++ (fasta s) ++ " -R "
           (rin, rout, rerr, rpid) <- runInteractiveCommand $ " reverse_complement.pl "

           extractStopCodon sinF sinR s annotationList

           f <- hGetContents soutF
           r <- hGetContents soutR
           hPutStr rin r
           reversed <- hGetContents rout
           putStr f
           putStr reversed
           waitForProcess pidF
           waitForProcess pidR
           waitForProcess rpid

           return ()

extractStopCodon:: Handle -> Handle -> Settings -> [Annotation] -> IO()
extractStopCodon forward reverse settings [] = return ()
extractStopCodon forward reverse settings (a:rest) = do
                                               saveStopCodonFromAnnotation forward reverse settings a
                                               extractStopCodon forward reverse settings rest
                                               return ()



saveStopCodonFromAnnotation :: Handle -> Handle -> Settings -> Annotation -> IO()
saveStopCodonFromAnnotation forward reverse s a = do
                                                     mapM (saveStopCodonFromGene forward reverse s) $ genes a
                                                     return ()

saveStopCodonFromGene :: Handle -> Handle -> Settings -> Gene -> IO()
saveStopCodonFromGene forward reverse s g = do
                                                mapM (saveStopCodonFromTranscript forward reverse s) $ transcripts g
                                                return ()


saveStopCodonFromTranscript :: Handle -> Handle -> Settings -> Transcript -> IO()
saveStopCodonFromTranscript forward reverse s t = do
                                      mapM (saveStopCodon forward reverse s t)  $ getStopCodonSites t
                                      return ()

saveStopCodon :: Handle -> Handle -> Settings -> Transcript -> Site -> IO()
saveStopCodon forward reverse s t site = case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s) + 1)
                                                             ((position site) - (offset s) + (length' s) )
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) )
                                                              ((position site) + (offset s) - 1)

printSequence ::  Handle -> Name -> Integer -> Integer -> IO()
printSequence  cdbyank n s e = do
                                 hPutStrLn  cdbyank  ( n ++ " " ++ (show s) ++ " " ++ (show e))
                                 return ()





