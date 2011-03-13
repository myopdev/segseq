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
import qualified Data.ByteString.Lazy as L
import Data.ByteString.Internal
import Directory
import FastaDB
import Data.List.Split

main :: IO()
main =  do
           s <-  (getArgs >>= (return . compilerOpts) >>= (liftM (buildSettings (Settings "" "" "" 9 3 (-1) (-1)))))
           gtfh <- openFile (gtf s) ReadMode
           indexFastaFile (fasta s) ((fasta s) ++ ".db")
           (B.hGetContents gtfh) >>= (return . createAnnotationList . readCSV . B.unpack) >>= (run s)
           hClose gtfh
           return ()

revcomp :: L.ByteString -> IO(L.ByteString)
revcomp seq = return (L.foldl complement L.empty (L.reverse seq))
 where complement y x            | (x == (c2w 'A')) = c2w 'T' `L.cons` y
                      | (x == (c2w 'a'))=  c2w 't' `L.cons` y
                      | (x == (c2w 'C'))=  c2w 'G' `L.cons` y
                      | (x == (c2w 'c')) = c2w 'g' `L.cons` y
                      | (x == (c2w 'G')) = c2w 'C' `L.cons` y
                      | (x == (c2w 'g')) = c2w 'c' `L.cons` y
                      | (x == (c2w 'T')) = c2w 'A' `L.cons` y
                      | (x == (c2w 't')) = c2w 'a' `L.cons` y
                      | otherwise   =   x `L.cons` y

printFeature :: Settings -> ([String],[String]) -> IO()
printFeature s (f,r)  =
  do forM f (\ seq -> getSequence (fasta s) ((fasta s) ++ ".db") seq >>= putSequence seq )
     forM r (\ seq -> getSequence (fasta s) ((fasta s) ++ ".db") seq >>=  revcomp >>= putSequence seq )
     return()
 where
putSequence :: String -> L.ByteString -> IO()
putSequence seq bytes =
  do if (bytes /= L.empty)
        then do putStrLn $ ">" ++ head (splitOn " " seq)
                L.putStrLn bytes
        else do return()



run :: Settings -> [Annotation] -> IO ()
run s a=   do
             features <- return(extractContent s a)
             printFeature s features
             return ();
