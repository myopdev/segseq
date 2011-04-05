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
import qualified GTF as GTF
import qualified GFF as GFF
import Data.Maybe
import System.Directory(getTemporaryDirectory, removeFile)
import Control.Concurrent
import qualified Data.ByteString.Char8 as B
import Directory
main :: IO()
main =  do
           args <- getArgs 
           gffh <- openFile (args !! 0) ReadMode
           (B.hGetContents gffh) >>= (return . GFF.createAnnotationList . GFF.readGFF . B.unpack) >>= (run)
           hClose gffh
           return ()


run ::  [Annotation] -> IO ()
run a =
  do forM a (putStrLn . GTF.renderGTF)
     return ();
