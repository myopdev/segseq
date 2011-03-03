module SegSequence where
import System.Cmd
import System.Console.GetOpt
import Data.Maybe (fromMaybe,fromJust)


data Settings = Settings  { gtf :: String,
                            fasta :: String,
                            length' :: Integer,
                            offset :: Integer }


data Flag = Version
          | GTF String
          | FASTA String
          | Output String
          | Length String
          | Offset String
          deriving Show


options :: [OptDescr Flag]
options =
  [  Option ['g'] ["gtf"] (ReqArg GTF "FILE") "gtf file name"
  ,  Option ['f'] ["fasta"] (ReqArg FASTA "FILE") "fasta file name"
  ,  Option ['l'] ["length"] (OptArg lengthp "Int") "length of the window"
  ,  Option ['o'] ["offset"] (OptArg offsetp "Int") "offset of the window" ]

lengthp :: Maybe String -> Flag
lengthp = Length . fromMaybe "9"
offsetp = Offset . fromMaybe "3"



compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
  case getOpt Permute options argv of
     (o,n,[]) -> return (o,n)
     (_,_,errs) -> ioError(userError (concat errs ++ usageInfo header options))
  where header = "Usage: segseq -g <gtf file> -f <fasta file> -d <output dir>"



createIndexFromFasta :: String  -> IO()
createIndexFromFasta s = do
                            io <- system  $ " cdbfasta " ++ s ++ " 2> cdbfasta.err 1> cdbfasta.out "
                            return ()



buildSettings :: Settings -> ([Flag],[String]) -> Settings
buildSettings settings (opts,n)  =  foldr nextOption defaults opts
       where nextOption option settings = case option of
                                                 GTF x -> Settings x
                                                                   (fasta settings)
                                                                   (length' settings)
                                                                   (offset settings)
                                                 FASTA y -> Settings (gtf settings)
                                                                     y
                                                                     (length' settings)
                                                                     (offset settings)
                                                 Length l -> Settings (gtf settings)
                                                                      (fasta settings)
                                                                      ((read (n!!0))  ::Integer)
                                                                      (offset settings)
                                                 Offset o -> Settings (gtf settings)
                                                                      (fasta settings)
                                                                      (length' settings)
                                                                      ((read  (n!!1)) :: Integer)
             defaults = (Settings "" ""  9 3)

