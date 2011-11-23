module Main where

import System
import SegSeq.FastaDB
import Data.Char
import Data.List
import System.IO

main :: IO()
main =  do
           s <-  getArgs
           indexFastaFile (head s) ((head s) ++ ".db")
           return ()

