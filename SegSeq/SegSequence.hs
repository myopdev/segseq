module SegSeq.SegSequence where
import System.Cmd
import System.Console.GetOpt
import Data.Maybe (fromMaybe,fromJust)
import Control.Monad
import System.Directory(getTemporaryDirectory, removeFile)
import System.IO
import SegSeq.SequenceAnnotation
import Data.List.Split

data Settings = Settings  { gtf :: String,
                            fasta :: String,
                            feature ::String,
                            length' :: Integer,
                            offset :: Integer,
                            phase' :: Integer,
                            size' :: Integer,
                            gc1 :: Integer,
                            gc2 :: Integer,
                            genefilter :: Bool,
                            clean::Bool,
                            alternative:: Bool  }

data Flag = Version
          | GTF String
          | FASTA String
          | Output String
          | Feature String
          | Length String
          | Offset String
          | Phase String
          | Size String
          | GC String
          | Remove String
          | RemoveAlternative String
          deriving Show

options :: [OptDescr Flag]
options =
  [  Option ['g'] ["gtf"] (ReqArg GTF "FILE") "gtf file name"
  ,  Option ['f'] ["fasta"] (ReqArg FASTA "FILE") "fasta file name"
  ,  Option ['x'] ["feature"] (ReqArg Feature "STRING") "the feature to extract"
  ,  Option ['l'] ["length"] (OptArg lengthp "Int") "length of the window"
  ,  Option ['o'] ["offset"] (OptArg offsetp "Int") "offset of the window"
  ,  Option ['p'] ["phase"] (OptArg phasep "Int") "get features in a specific phase"
  ,  Option ['s'] ["size"] (OptArg sizep "Int") "get features with size at most s"
  ,  Option ['c'] ["gc content"] (OptArg gcp "STRING") "filter by gc"
  ,  Option ['r'] ["remove invalids"] (OptArg qp "Bool") "remove sequences that contain strange nucleotides symbols and exclude duplicates "
  ,  Option ['a'] ["only one variant"] (OptArg alternativep "Bool") "get features from only one transcript variant per gene " ]

lengthp,offsetp,gcp,qp,alternativep :: Maybe String -> Flag
lengthp = Length . fromMaybe "9"
offsetp = Offset . fromMaybe "3"
phasep = Phase . fromMaybe "-1"
sizep = Size . fromMaybe "-1"
gcp = GC . fromMaybe "-1:-1"
qp = Remove . fromMaybe "-1"
alternativep = RemoveAlternative . fromMaybe "-1"


compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
  case getOpt Permute options argv of
     (o,n,errs) -> case (length o == 0) of
                     True -> ioError(userError (concat errs ++ usageInfo header options))
                     False -> case (length errs == 0) of
                                True -> return (o,n)
                                False-> ioError(userError (concat errs ++ usageInfo header options))
  where header = "Usage: segseq -g <gtf file> -f <fasta file> -x [initial|internal|final|all-exons|intron|sites] [-l <length of the window>] [-o <offset>] [-s <size>]"

buildSettings :: Settings -> ([Flag],[String]) -> Settings
buildSettings settings (opts,n)  =  fst (foldl nextOption (defaults,0) opts)
       where nextOption (Settings g fa fe le o p s g1 g2 f q var,count)  option  =
                                         case option of
                                                 GTF x -> (Settings x fa fe le o p s g1 g2 f q var, count)
                                                 FASTA x -> (Settings g x fe le o p s  g1 g2 f q var, count)
                                                 Feature x ->(Settings g fa x le o p s g1 g2 f q var, count)
                                                 Length x ->  (Settings g fa fe ((read (n!!count))  ::Integer) o p s g1 g2 f q var, count + 1)
                                                 Offset x ->(Settings g fa fe le ((read (n!!count))  ::Integer) p s  g1 g2 f q var, count + 1)
                                                 Phase x -> (Settings g fa fe le o ((read (n!!count))  ::Integer) s  g1 g2 f q var, count + 1)
                                                 Size x -> (Settings g fa fe le o p ((read (n!!count))  ::Integer)  g1 g2 f q var, count + 1)
                                                 GC x -> (Settings g fa fe le o p s  (getg1 (n!!count)) (getg2 (n!!count)) (getrestriction (n!!count)) q var, count + 1)
                                                 Remove x ->((Settings g fa fe le o p s g1 g2 f True var), (count))
                                                 RemoveAlternative x ->((Settings g fa fe le o p s g1 g2 f q True), (count))
             defaults = Settings "" "" ""  9 3 (-1) (-1) (-1) (-1) False False False
             getg1 s = read ((splitOn ":" s)!!0) :: Integer
             getg2 s = read ((splitOn ":" s)!!1) :: Integer
             getrestriction s = case (length (splitOn ":" s) >= 3) of
                                   True -> ((splitOn ":" s) !!2) == "gene"
                                   False -> False


extractFeature :: Settings -> [a] -> (Settings -> a -> [b] -> ([String],[String]) ) -> (a ->[b]) -> ([String], [String])
extractFeature s [] k y = ([],[])
extractFeature s (a:rest) k y = (fst(r1) ++ fst(r2), snd(r1) ++ snd(r2))
                               where  r1 = k s a (y a)
                                      r2 = extractFeature s rest k y


extractFeatureList :: Settings -> [a] -> (Settings -> a -> [b] -> [([String],[String])] ) -> (a ->[b]) -> [([String], [String])]
extractFeatureList s [] k y = [([],[])]
extractFeatureList s (a:rest) k y = r1 ++ r2
                               where  r1 = k s a (y a)
                                      r2 = extractFeatureList s rest k y


extractFeatureJoin :: Settings -> [a] -> (Settings -> a -> [b] -> ([String],[String]) ) -> (a ->[b]) -> ([String], [String])
extractFeatureJoin s [] k y = ([],[])
extractFeatureJoin s (a:rest) k y =  ([(head' f1) ++ (head' f2)], [(head' r1) ++ (head' r2)])
                               where  (f1,r1) = k s a (y a)
                                      (f2,r2) = extractFeatureJoin s rest k y


head' :: [String] -> String
head' l = case (length l >= 1) of
               True -> head l
               False -> ""
extractContent :: Settings -> [Annotation] -> ([String], [String])
extractContent s a | feature s == "initial" = extractInitialExons s a
                   | feature s == "internal" = extractInternalExons s a
                   | feature s == "final" = extractFinalExons s a
                   | feature s == "acceptor" = extractSites s a
                   | feature s == "donor" = extractSites s a
                   | feature s == "start" = extractSites s a
                   | feature s == "stop" = extractSites s a
                   | feature s == "intron" = extractIntrons s a
                   | feature s == "exons" = extractExons s a
                   | feature s == "intergenic" = extractIntergenic s a
                   | feature s == "single" = extractSingleExons s a
                   | feature s == "cds" = extractCDS s a
                   | feature s == "noncoding" = extractNonCoding s a
                   | feature s == "initial-pattern" = extractInitialPattern s a
                   | otherwise = extractExons s a




extractIntergenic :: Settings -> [Annotation] -> ([String], [String])
extractIntergenic s a = (concat $ map (\ annot -> printIntergenic annot) a, [])
                  where printIntergenic annot = map (\ positions -> printSequence s name (fst positions) (snd positions) ) (getIntergenicRegions annot)
                          where name = seqname (seqentry annot)


extractNonCoding :: Settings -> [Annotation] -> ([String], [String])
extractNonCoding s a = (concat $ map (\ annot -> printIntergenic annot) a, [])
                  where printIntergenic annot = map (\ positions -> printSequence s name (fst positions) (snd positions) ) (getNonCodingRegions annot)
                          where name = seqname (seqentry annot)



extractIntrons :: Settings -> [Annotation] -> ([String], [String])
extractIntrons s a = extractFeature s a fromGenes (\ annot -> genes annot)
                    where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                          fromTranscripts s p t = extractFeature s t fromIntrons ( \ txs -> (getIntrons txs) )
                          fromIntrons s p introns = case (strand p == "+") of
                                                       True ->  (seqstr  , [])
                                                       False -> ([],seqstr)
                                                    where seqstr = map ( \ i -> (printSequence s seqname (fst i) (snd i)) ++ " " ++ (txname p)) introns
                                                          seqname = parentseq (cstart ((txcds p) !! 0))




extractSingleExons :: Settings -> [Annotation] -> ([String], [String])
extractSingleExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                       fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> txcds x ) 
                       fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                       fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> case (stype (cend c)) of
                                                            StopCodon -> forward
                                                            _ -> ([],[])
                                             StopCodon -> case (stype (cend c)) of
                                                           StartCodon -> reverse
                                                           _ -> ([],[])
                                             _ -> ([],[])
                                             where forward = ([seqstr],[])
                                                   reverse = ([], [seqstr])
                                                   seqstr = (printSequence s (parentseq (cstart c))
                                                                          (position (cstart c))
                                                                          (position (cend c))  ++ " " ++ (parentTxname (cstart c)))



extractFinalExons :: Settings -> [Annotation] -> ([String], [String])
extractFinalExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                       fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithStartPhase (phase' s) x) )
                       fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                       fromCDS s p (c:_)= case (stype (cstart c)) of
                                             Acceptor -> case (stype (cend c)) of
                                                           StopCodon -> forward
                                                           _ -> ([],[])
                                             StartCodon -> case (stype (cend c)) of
                                                           StopCodon -> forward
                                                           _ -> ([],[])
                                             StopCodon -> case (stype (cend c)) of
                                                           Acceptor -> reverse
                                                           _ -> ([],[])
                                             _ -> ([],[])
                                             where forward = ([seqstr],[])
                                                   reverse = ([], [seqstr])
                                                   seqstr = (printSequence s (parentseq (cstart c))
                                                                          (position (cstart c))
                                                                          (position (cend c))  ++ " " ++ (parentTxname (cstart c)))



extractInternalExons :: Settings -> [Annotation] -> ([String], [String])
extractInternalExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                       fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithStartPhase (phase' s) x))
                       fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                       fromCDS s p (c:_)= case (stype (cstart c)) of
                                             Acceptor -> case (stype (cend c)) of
                                                           Donor -> forward
                                                           _ -> ([],[])
                                             Donor -> case (stype (cend c)) of
                                                           Acceptor -> reverse
                                                           _ -> ([],[])
                                             _ -> ([],[])
                                             where forward = ([seqstr],[])
                                                   reverse = ([], [seqstr])
                                                   seqstr = (printSequence s (parentseq (cstart c))
                                                                          (position (cstart c))
                                                                          (position (cend c)) ++ " " ++ (parentTxname (cstart c)))



extractInitialExons :: Settings -> [Annotation] -> ([String], [String])
extractInitialExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                       fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithEndPhase (phase' s) x))
                       fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                       fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> forward
                                             Donor -> case (stype (cend c)) of
                                                        StartCodon -> reverse
                                                        _ -> ([],[])
                                             StopCodon -> case (stype (cend c)) of
                                                            StartCodon -> reverse
                                                            _ -> ([],[])

                                             _ -> ([],[])
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstr])
                                                  seqstr = (printSequence  s (parentseq (cstart c))
                                                                         (position (cstart c))
                                                                         (position (cend c)) ++ " " ++ (parentTxname (cstart c)))




extractCDS :: Settings -> [Annotation] -> ([String], [String])
extractCDS s a = extractFeature s a fromGenes (\ x -> genes x)
                   where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                         fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> txcds x )
                         fromCDSList s p c = ([case (length (fwd')>= 1 && (head fwd') /= "") of
                                                    True -> (parentseq (cstart (c!!0))) ++ " " ++ (head fwd')++"\n"
                                                    False -> ""],
                                              [case (length (rev') >= 1 && (head rev') /= "") of
                                                    True -> (parentseq (cstart (c!!0))) ++ " " ++ (head rev') ++ "\n"
                                                    False -> ""])
                                            where (fwd',rev') = extractFeatureJoin s c fromCDS ( \x -> [x] )

                         fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> forward
                                             Acceptor -> forward
                                             Donor -> reverse
                                             StopCodon -> reverse
                                             Unknown -> ([],[])
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstr])
                                                  seqstr = " "++ show (position (cstart c)) ++ " " ++ show  (position (cend c))



extractExons :: Settings -> [Annotation] -> ([String], [String])
extractExons s a = extractFeature s a fromGenes (\ x -> genes x)
                   where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                         fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithStartPhase (phase' s) x))
                         fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                         fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> forward
                                             Acceptor -> forward
                                             Donor -> reverse
                                             StopCodon -> reverse
                                             Unknown -> ([], [])
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstr])
                                                  seqstr = (printSequence s  (parentseq (cstart c))
                                                                         (position (cstart c))
                                                                         (position (cend c)) ++ " " ++ (parentTxname (cstart c)))




extractInitialPattern :: Settings -> [Annotation] -> ([String], [String])
extractInitialPattern s a = extractFeature s a fromGenes (\ x -> genes x)
                   where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                         fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithStartPhase (phase' s) x))
                         fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                         fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> forward
                                             Acceptor -> forward
                                             Donor -> reverse
                                             StopCodon -> reverse
                                             Unknown -> ([], [])
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstrR])
                                                  seqstr = (printSequence s  (parentseq (cstart c))
                                                                            (position (cstart c) + (offset s))
                                                                            (position (cstart c) + (offset s) + (length' s) - 1) ++ " " ++ (parentTxname (cstart c)))


                                                  seqstrR = (printSequence s (parentseq (cend c))
                                                                            (position (cend c)  - (offset s))
                                                                            (position (cend c) - (offset s) - (length' s) + 1)++ " " ++ (parentTxname (cstart c)))



extractSites :: Settings -> [Annotation] -> ([String], [String])
extractSites s a = extractFeature s a fromGenes (\ annot -> genes annot)
                    where fromGenes s p g = extractFeature s g fromTranscripts (filterTranscript s)
                          fromTranscripts s p t = extractFeature s t fromSiteList ( \ txs -> (getSites (phase' s) (siteName s) txs) )
                          fromSiteList s p sites = getSiteString s p sites


getSiteString ::  Settings -> Transcript -> [Site] -> ([String], [String])
getSiteString s t [] =  ([],[] )
getSiteString s t (site:sites) = (fst (x) ++ fst(y), snd(x) ++ snd(y))
                                  where x = singleSiteToStr s t site
                                        y = getSiteString s t sites

singleSiteToStr :: Settings -> Transcript -> Site -> ([String], [String])
singleSiteToStr s t site =
         case stype site of
            StartCodon ->case (strand t == "+") of
                               True -> ([(printSequence  s (parentseq site)
                                                        ((position site) - (offset s))
                                                        ((position site) - (offset s) + (length' s) - 1)) ++ " " ++ (parentTxname site)],[])
                               False -> ([],[printSequence s (parentseq site)
                                                           ((position site) + (offset s) - (length' s) + 1)
                                                           ((position site) + (offset s)) ++ " " ++ (parentTxname  site)])
            Donor -> case (strand t == "+") of
                       True -> ([printSequence s (parentseq site)
                                             ((position site) - (offset s) + 1)
                                             ((position site) - (offset s) + (length' s) )++ " " ++ (parentTxname  site)] ,[])
                       False -> ([], [printSequence s (parentseq site)
                                                    ((position site) + (offset s) - (length' s) )
                                              ((position site) + (offset s) - 1)++ " " ++ (parentTxname  site)])

            Acceptor ->  case (strand t == "+") of
                               True -> ([printSequence s (parentseq site)
                                                       ((position site) - (offset s) - 2)
                                                       ((position site) - (offset s) + (length' s) - 3)++ " " ++ (parentTxname  site)],[])
                               False -> ([], [printSequence s  (parentseq site)
                                                            ((position site) + (offset s) - (length' s) + 3)
                                                             ((position site) + (offset s) + 2) ++ " " ++ (parentTxname  site)])

            StopCodon -> case (strand t == "+") of
                               True -> ([printSequence s (parentseq site)
                                                       ((position site) - (offset s) + 1)
                                                       ((position site) - (offset s) + (length' s) )++ " " ++ (parentTxname  site)], [])
                               False -> ([], [printSequence s (parentseq site)
                                                            ((position site) + (offset s) - (length' s) )
                                                            ((position site) + (offset s) - 1)++ " " ++ (parentTxname  site)])

siteName :: Settings -> SiteType
siteName s | feature s == "start" = StartCodon
           | feature s == "stop" = StopCodon
           | feature s == "acceptor" = Acceptor
           | feature s == "donor" = Donor
           | otherwise = StartCodon



printSequence :: Settings -> Name -> Integer -> Integer -> String
printSequence set n s e = case ((size' set) > 0) of
                               True -> case  ((e -s + 1) <= (size' set)) of
                                             True -> n ++ " " ++ show (s) ++ " " ++ show e ++ "\n"
                                             False -> ""
                               False -> n ++ " " ++ show (s) ++ " " ++ show e ++ "\n"


filterTranscript :: Settings -> Gene -> [Transcript]
filterTranscript s g = case (alternative s) of
                         True -> [(transcripts g )!!0]
                         False -> transcripts g
