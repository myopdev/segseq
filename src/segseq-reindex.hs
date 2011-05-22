module Main where

import System
import FastaDB
import Data.Char
import Data.List
import FastaDB
import System.IO

main :: IO()
main =  do
           s <-  getArgs
           indexFastaFile (head s) ((head s) ++ ".db")
           return ()

