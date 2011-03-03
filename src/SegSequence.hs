module SegSequence where
import System.Cmd
import System.Console.GetOpt
import Data.Maybe (fromMaybe,fromJust)
import Control.Monad
import System.Directory(getTemporaryDirectory, removeFile)
import System.IO
import SequenceAnnotation

data Settings = Settings  { gtf :: String,
                            fasta :: String,
                            feature ::String,
                            length' :: Integer,
                            offset :: Integer }


data Flag = Version
          | GTF String
          | FASTA String
          | Output String
          | Feature String
          | Length String
          | Offset String
          deriving Show

options :: [OptDescr Flag]
options =
  [  Option ['g'] ["gtf"] (ReqArg GTF "FILE") "gtf file name"
  ,  Option ['f'] ["fasta"] (ReqArg FASTA "FILE") "fasta file name"
  ,  Option ['x'] ["feature"] (ReqArg Feature "STRING") "the feature to extract"
  ,  Option ['l'] ["length"] (OptArg lengthp "Int") "length of the window"
  ,  Option ['o'] ["offset"] (OptArg offsetp "Int") "offset of the window" ]

lengthp,offsetp :: Maybe String -> Flag
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
                                                                   (feature settings)
                                                                   (length' settings)
                                                                   (offset settings)
                                                 FASTA y -> Settings (gtf settings)
                                                                     y
                                                                     (feature settings)
                                                                     (length' settings)
                                                                     (offset settings)
                                                 Length l -> Settings (gtf settings)
                                                                      (fasta settings)
                                                                      (feature settings)
                                                                      ((read (n!!0))  ::Integer)
                                                                      (offset settings)
                                                 Offset o -> Settings (gtf settings)
                                                                      (fasta settings)
                                                                      (feature settings)
                                                                      (length' settings)
                                                                      ((read  (n!!1)) :: Integer)
                                                 Feature f -> Settings (gtf settings)
                                                                       (fasta settings)
                                                                       f
                                                                       (length' settings)
                                                                       (offset settings)

             defaults = (Settings "" "" ""  9 3)





extractSite:: SiteType->Handle -> Handle -> Settings -> [Annotation] -> IO()
extractSite stype' forward reverse settings  [] = return ()
extractSite stype' forward reverse settings  (a:rest) = do
                                                         saveSiteFromAnnotation stype' forward reverse settings a
                                                         extractSite stype' forward reverse settings rest
                                                         return ()



saveSiteFromAnnotation :: SiteType -> Handle -> Handle -> Settings ->  Annotation -> IO()
saveSiteFromAnnotation stype' forward reverse s a  = do  mapM (saveSiteFromGene stype' forward reverse s) $ genes a
                                                         return ()

saveSiteFromGene :: SiteType -> Handle -> Handle -> Settings -> Gene -> IO()
saveSiteFromGene stype' forward reverse s g = do
                                                mapM (saveSiteFromTranscript stype' forward reverse s) $ transcripts g
                                                return ()
saveSiteFromTranscript :: SiteType -> Handle -> Handle -> Settings -> Transcript -> IO()
saveSiteFromTranscript stype' forward reverse s  t = do
                                                      mapM (saveSite forward reverse s t)  $ (getSites stype') t
                                                      return ()


printSequence ::  Handle -> Name -> Integer -> Integer -> IO()
printSequence  cdbyank n s e = do
                                 hPutStrLn  cdbyank  ( n ++ " " ++ (show s) ++ " " ++ (show e))
                                 return ()


saveSite :: Handle -> Handle -> Settings -> Transcript -> Site -> IO()
saveSite forward reverse s t site =
         case stype site of
            StartCodon ->case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s))
                                                             ((position site) - (offset s) + (length' s) - 1)
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) + 1)
                                                              ((position site) + (offset s))
            Donor -> case (strand t == "+") of
                       True -> printSequence forward (parentseq site)
                                                     ((position site) - (offset s) + 1)
                                                     ((position site) - (offset s) + (length' s) )
                       False -> printSequence reverse (parentseq site)
                                                      ((position site) + (offset s) - (length' s) )
                                                      ((position site) + (offset s) - 1)

            Acceptor ->  case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s) - 2)
                                                             ((position site) - (offset s) + (length' s) - 3)
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) + 3)
                                                              ((position site) + (offset s) + 2)

            StopCodon -> case (strand t == "+") of
                               True -> printSequence forward (parentseq site)
                                                             ((position site) - (offset s) + 1)
                                                             ((position site) - (offset s) + (length' s) )
                               False -> printSequence reverse (parentseq site)
                                                              ((position site) + (offset s) - (length' s) )
                                                              ((position site) + (offset s) - 1)

siteName :: Settings -> SiteType
siteName s | feature s == "start" = StartCodon
           | feature s == "stop" = StopCodon
           | feature s == "acceptor" = Acceptor
           | feature s == "donor" = Donor
           | otherwise = StartCodon