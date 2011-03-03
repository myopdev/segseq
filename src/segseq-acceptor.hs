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

           extractAcceptorSite sinF sinR s annotationList

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


extractAcceptorSite:: Handle -> Handle -> Settings -> [Annotation] -> IO()
extractAcceptorSite forward reverse settings [] = return ()
extractAcceptorSite forward reverse settings (a:rest) = do
                                               saveAcceptorSiteFromAnnotation forward reverse settings a
                                               extractAcceptorSite forward reverse settings rest
                                               return ()



saveAcceptorSiteFromAnnotation :: Handle -> Handle -> Settings -> Annotation -> IO()
saveAcceptorSiteFromAnnotation forward reverse s a = do  mapM (saveAcceptorSiteFromGene forward reverse s) $ genes a
                                                         return ()

saveAcceptorSiteFromGene :: Handle -> Handle -> Settings -> Gene -> IO()
saveAcceptorSiteFromGene forward reverse s g = do
                                                mapM (saveAcceptorSiteFromTranscript forward reverse s) $ transcripts g
                                                return ()


saveAcceptorSiteFromTranscript :: Handle -> Handle -> Settings -> Transcript -> IO()
saveAcceptorSiteFromTranscript forward reverse s t = do
                                                      mapM (saveAcceptorSite forward reverse s t)  $ getAcceptorSites t
                                                      return ()

saveAcceptorSite :: Handle -> Handle -> Settings -> Transcript -> Site -> IO()
saveAcceptorSite forward reverse s t site = case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s) - 2)
                                                             ((position site) - (offset s) + (length' s) - 3)
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) + 3)
                                                              ((position site) + (offset s) + 2)

printSequence ::  Handle -> Name -> Integer -> Integer -> IO()
printSequence  cdbyank n s e = do
                                 hPutStrLn  cdbyank  ( n ++ " " ++ (show s) ++ " " ++ (show e))
                                 return ()



