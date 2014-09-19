module Main where

import SegSeq.FastaDB
import Data.Char
import Data.List
import System.IO
import System.Environment (getArgs)

main :: IO()
main =  do
           s <-  getArgs
           indexFastaFile (head s) ((head s) ++ ".db")
           return ()

