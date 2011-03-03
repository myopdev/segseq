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

           extractDonorSite sinF sinR s annotationList

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


extractDonorSite:: Handle -> Handle -> Settings -> [Annotation] -> IO()
extractDonorSite forward reverse settings [] = return ()
extractDonorSite forward reverse settings (a:rest) = do
                                               saveDonorSiteFromAnnotation forward reverse settings a
                                               extractDonorSite forward reverse settings rest
                                               return ()



saveDonorSiteFromAnnotation :: Handle -> Handle -> Settings -> Annotation -> IO()
saveDonorSiteFromAnnotation forward reverse s a = do  mapM (saveDonorSiteFromGene forward reverse s) $ genes a
                                                      return ()

saveDonorSiteFromGene :: Handle -> Handle -> Settings -> Gene -> IO()
saveDonorSiteFromGene forward reverse s g = do
                                                mapM (saveDonorSiteFromTranscript forward reverse s) $ transcripts g
                                                return ()


saveDonorSiteFromTranscript :: Handle -> Handle -> Settings -> Transcript -> IO()
saveDonorSiteFromTranscript forward reverse s t = do
                                                      mapM (saveDonorSite forward reverse s t)  $ getDonorSites t
                                                      return ()

saveDonorSite :: Handle -> Handle -> Settings -> Transcript -> Site -> IO()
saveDonorSite forward reverse s t site = case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s)+1)
                                                             ((position site) - (offset s) + (length' s) )
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) )
                                                              ((position site) + (offset s) -1)

printSequence ::  Handle -> Name -> Integer -> Integer -> IO()
printSequence  cdbyank n s e = do
                                 hPutStrLn  cdbyank  ( n ++ " " ++ (show s) ++ " " ++ (show e))
                                 return ()



