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
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Lazy as L
import Data.ByteString.Internal
import System.Directory
import FastaDB
import Data.List.Split
import Data.Char
import Data.List

main :: IO()
main =  do
           s <-  (getArgs >>= (return . compilerOpts) >>= (liftM (buildSettings (Settings "" "" "" 9 3 (-1) (-1) (-1) (-1) (False) (False) (False)))))
           hPutStrLn stderr $ "r: " ++ (show $ clean s)
           hPutStrLn stderr $ "x: " ++ (feature s)
           gtfh <- openFile (gtf s) ReadMode
           indexFastaFile (fasta s) ((fasta s) ++ ".db")
           (B.hGetContents gtfh) >>= (return . createAnnotationList . readCSV . B.unpack) >>= (filterGeneByGC s) >>= (run s)
           hClose gtfh
           return ()

revcomp :: L.ByteString -> IO(L.ByteString)
revcomp seq = return (L.foldl' complement L.empty (seq))
 where complement y x | (x == (c2w 'A')) = c2w 'T' `L.cons` y
                      | (x == (c2w 'a'))=  c2w 't' `L.cons` y
                      | (x == (c2w 'C'))=  c2w 'G' `L.cons` y
                      | (x == (c2w 'c')) = c2w 'g' `L.cons` y
                      | (x == (c2w 'G')) = c2w 'C' `L.cons` y
                      | (x == (c2w 'g')) = c2w 'c' `L.cons` y
                      | (x == (c2w 'T')) = c2w 'A' `L.cons` y
                      | (x == (c2w 't')) = c2w 'a' `L.cons` y
                      | otherwise   =   x `L.cons` y

data SeqComposition = SeqComposition !Int !Int !Int !Int

getComposition :: L.ByteString -> IO(SeqComposition)
getComposition seq = return $! (L.foldl' sumComposition (SeqComposition 0 0 0 0)  (seq))
 where sumComposition (SeqComposition a c g t)  x | (x == (c2w 'A')) = SeqComposition (a+1) c g t
                                                  | (x == (c2w 'a'))=  SeqComposition (a+1) c g t
                                                  | (x == (c2w 'C'))=  SeqComposition a (c+1) g t
                                                  | (x == (c2w 'c')) = SeqComposition a (c+1) g t
                                                  | (x == (c2w 'G')) = SeqComposition a c (g+1) t
                                                  | (x == (c2w 'g')) = SeqComposition a c (g+1) t
                                                  | (x == (c2w 'T')) = SeqComposition a c g (t+1)
                                                  | (x == (c2w 't')) = SeqComposition a c g (t+1)
                                                  | otherwise   =   SeqComposition a c g t



cleanDataset :: Settings -> L.ByteString -> IO(L.ByteString)
cleanDataset s  bytes =
  do if ( clean s )
        then do case (L.find invalid bytes) of
                     Nothing -> return(bytes)
                     Just n -> return(L.empty)
        else do return (bytes);
  where invalid w | ( w == (c2w 'A') ) = False
                  | ( w == (c2w 'a') ) = False
                  | ( w == (c2w 'C') ) = False
                  | ( w == (c2w 'c') ) = False
                  | ( w == (c2w 'G') ) = False
                  | ( w == (c2w 'g') ) = False
                  | ( w == (c2w 'T') ) = False
                  | ( w == (c2w 't') ) = False
                  | otherwise = True


filterByGC :: Settings -> L.ByteString -> IO(L.ByteString)
filterByGC s  bytes =
  do if ( ((gc1 s) >= 0) && ((gc2 s) >= 0))
        then do composition <- getComposition bytes
                gc <- return( getGC composition)
                if (((gc1 s) <= gc) && (gc <= (gc2 s)))
                   then do return (bytes)
                   else do return (L.empty)
        else do return(bytes)

getGC :: SeqComposition -> Integer
getGC (SeqComposition a c g t) = floor $ 100.0 * (fromIntegral (g + c)/ fromIntegral (a+c+g+t))



printFeature :: Settings -> ([String],[String]) -> IO()
printFeature s (f,r)  =
  do forM cleanF (\ seq -> getSequence (fasta s) ((fasta s) ++ ".db") seq >>=  filterByGC s  >>= cleanDataset s >>= putSequence seq )
     forM cleanR (\ seq -> getSequence (fasta s) ((fasta s) ++ ".db") seq >>=  revcomp >>=  filterByGC s >>= cleanDataset s >>= putSequence seq )
     return()
 where cleanF = case (clean s) of
                  True -> nub f
                  False -> f
       cleanR = case (clean s) of
                  True -> nub r
                  False -> r




putSequence :: String -> L.ByteString -> IO()
putSequence  seq bytes =
  do if (bytes /= L.empty)
        then do putStrLn $ ">" ++ head (splitOn " " seq)
                L.putStrLn bytes
        else do return()



run :: Settings -> [Annotation] -> IO ()
run s a=
 do features <- return(extractContent s a)
    printFeature s features
    return ();




filterGeneByGC :: Settings -> [Annotation] -> IO([Annotation])
filterGeneByGC s a =
  do if(genefilter s)
       then do foldM (buildNewAnnotationList s) [] a
       else do return (a)
 where buildNewAnnotationList s list (Annotation seqentry g source)  =
        do filtered <- foldM (filterGenes s (seqname $ seqentry)) [] g
           return (list ++ [Annotation seqentry filtered source])
       filterGenes s name list g =
        do positions <- return (getGenePosition g)
           seqEntry <- return (printSequence s name (fst positions) (snd positions))
           seq <- getSequence (fasta s) ((fasta s) ++ ".db") seqEntry >>= filterByGC s
           if (seq  /= L.empty)
              then do return([g] ++ list)
              else do return(list)







