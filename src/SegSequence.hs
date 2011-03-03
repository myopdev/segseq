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




extractSites:: SiteType -> Settings -> [Annotation] -> ([String], [String])
extractSites stype' settings  []  =  ([],[] )
extractSites  stype' settings  (a:rest)  =   (fst (x) ++ fst(y), snd(x) ++ snd(y))
                                             where x = sitesFromGenes stype' settings (genes a)
                                                   y = extractSites stype' settings rest

sitesFromGenes :: SiteType -> Settings -> [Gene] -> ([String], [String])
sitesFromGenes stype' s [] =  ([],[] )
sitesFromGenes stype' s (g:genes)  = (fst (x) ++ fst(y), snd(x) ++ snd(y))
                                      where x = sitesFromTranscripts stype' s (transcripts g)
                                            y = sitesFromGenes stype' s genes

sitesFromTranscripts :: SiteType -> Settings -> [Transcript] -> ([String], [String])
sitesFromTranscripts stype' s []  =   ([],[] )
sitesFromTranscripts stype' s (t:txs) = (fst (x) ++ fst(y), snd(x) ++ snd(y))
                                       where x = getSiteString s t (getSites stype' t)
                                             y = sitesFromTranscripts stype' s txs


getSiteString ::  Settings -> Transcript -> [Site] -> ([String], [String])
getSiteString s t [] =  ([],[] )
getSiteString  s t (site:sites) = (fst (x) ++ fst(y), snd(x) ++ snd(y))
                                  where x = singleSiteToStr s t site
                                        y = getSiteString s t sites




singleSiteToStr :: Settings -> Transcript -> Site -> ([String], [String])
singleSiteToStr s t site =
         case stype site of
            StartCodon ->case (strand t == "+") of
                               True -> ([printSequence  (parentseq site)
                                                        ((position site) - (offset s))
                                                        ((position site) - (offset s) + (length' s) - 1)],[])
                               False -> ([],[printSequence (parentseq site)
                                                           ((position site) + (offset s) - (length' s) + 1)
                                                           ((position site) + (offset s))])
            Donor -> case (strand t == "+") of
                       True -> ([printSequence (parentseq site)
                                             ((position site) - (offset s) + 1)
                                             ((position site) - (offset s) + (length' s) )] ,[])
                       False -> ([], [printSequence (parentseq site)
                                                    ((position site) + (offset s) - (length' s) )
                                              ((position site) + (offset s) - 1)])

            Acceptor ->  case (strand t == "+") of
                               True -> ([printSequence (parentseq site)
                                                       ((position site) - (offset s) - 2)
                                                       ((position site) - (offset s) + (length' s) - 3)],[])
                               False -> ([], [printSequence (parentseq site)
                                                            ((position site) + (offset s) - (length' s) + 3)
                                                             ((position site) + (offset s) + 2) ])

            StopCodon -> case (strand t == "+") of
                               True -> ([printSequence (parentseq site)
                                                       ((position site) - (offset s) + 1)
                                                       ((position site) - (offset s) + (length' s) )], [])
                               False -> ([], [printSequence (parentseq site)
                                                            ((position site) + (offset s) - (length' s) )
                                                            ((position site) + (offset s) - 1)])

siteName :: Settings -> SiteType
siteName s | feature s == "start" = StartCodon
           | feature s == "stop" = StopCodon
           | feature s == "acceptor" = Acceptor
           | feature s == "donor" = Donor
           | otherwise = StartCodon



printSequence :: Name -> Integer -> Integer -> String
printSequence n s e = n ++ " " ++ " " ++ show (s) ++ " " ++ show e