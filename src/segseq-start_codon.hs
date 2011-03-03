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

           extractStartCodon sinF sinR s annotationList

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


extractStartCodon:: Handle -> Handle -> Settings -> [Annotation] -> IO()
extractStartCodon forward reverse settings [] = return ()
extractStartCodon forward reverse settings (a:rest) = do
                                               saveStartCodonFromAnnotation forward reverse settings a
                                               extractStartCodon forward reverse settings rest
                                               return ()



saveStartCodonFromAnnotation :: Handle -> Handle -> Settings -> Annotation -> IO()
saveStartCodonFromAnnotation forward reverse s a = do  mapM (saveStartCodonFromGene forward reverse s) $ genes a
                                                       return ()

saveStartCodonFromGene :: Handle -> Handle -> Settings -> Gene -> IO()
saveStartCodonFromGene forward reverse s g = do
                                                mapM (saveStartCodonFromTranscript forward reverse s) $ transcripts g
                                                return ()


saveStartCodonFromTranscript :: Handle -> Handle -> Settings -> Transcript -> IO()
saveStartCodonFromTranscript forward reverse s t = do
                                                      mapM (saveStartCodon forward reverse s t)  $ getStartCodonSites t
                                                      return ()

saveStartCodon :: Handle -> Handle -> Settings -> Transcript -> Site -> IO()
saveStartCodon forward reverse s t site = case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s))
                                                             ((position site) - (offset s) + (length' s) - 1)
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) + 1)
                                                              ((position site) + (offset s))

printSequence ::  Handle -> Name -> Integer -> Integer -> IO()
printSequence  cdbyank n s e = do
                                 hPutStrLn  cdbyank  ( n ++ " " ++ (show s) ++ " " ++ (show e))
                                 return ()



