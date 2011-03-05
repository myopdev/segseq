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
                            offset :: Integer,
                            phase' :: Integer}

data Flag = Version
          | GTF String
          | FASTA String
          | Output String
          | Feature String
          | Length String
          | Offset String
          | Phase String
          deriving Show

options :: [OptDescr Flag]
options =
  [  Option ['g'] ["gtf"] (ReqArg GTF "FILE") "gtf file name"
  ,  Option ['f'] ["fasta"] (ReqArg FASTA "FILE") "fasta file name"
  ,  Option ['x'] ["feature"] (ReqArg Feature "STRING") "the feature to extract"
  ,  Option ['l'] ["length"] (OptArg lengthp "Int") "length of the window"
  ,  Option ['o'] ["offset"] (OptArg offsetp "Int") "offset of the window"
  ,  Option ['p'] ["phase"] (OptArg phasep "Int") "get features in a specific phase" ]

lengthp,offsetp :: Maybe String -> Flag
lengthp = Length . fromMaybe "9"
offsetp = Offset . fromMaybe "3"
phasep = Phase . fromMaybe "-1"

compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
  case getOpt Permute options argv of
     (o,n,errs) -> case (length o == 0) of
                     True -> ioError(userError (concat errs ++ usageInfo header options))
                     False -> case (length errs == 0) of
                                True -> return (o,n)
                                False-> ioError(userError (concat errs ++ usageInfo header options))
  where header = "Usage: segseq -g <gtf file> -f <fasta file> -x [initial|internal|final|all-exons|intron|sites]"

buildSettings :: Settings -> ([Flag],[String]) -> Settings
buildSettings settings (opts,n)  =  fst (foldl nextOption (defaults,0) opts)
       where nextOption (settings,count)  option  =
                                         case option of
                                                 GTF x -> (Settings x
                                                                   (fasta settings)
                                                                   (feature settings)
                                                                   (length' settings)
                                                                   (offset settings)
                                                                   (phase' settings), count)
                                                 FASTA y -> (Settings (gtf settings)
                                                                      y
                                                                      (feature settings)
                                                                      (length' settings)
                                                                      (offset settings)
                                                                      (phase' settings), count)
                                                 Length l -> (Settings (gtf settings)
                                                                       (fasta settings)
                                                                       (feature settings)
                                                                       ((read (n!!count))  ::Integer)
                                                                       (offset settings)
                                                                       (phase' settings), count + 1)
                                                 Offset o -> (Settings (gtf settings)
                                                                       (fasta settings)
                                                                       (feature settings)
                                                                       (length' settings)
                                                                       ((read  (n!!count)) :: Integer)
                                                                       (phase' settings), count + 1)
                                                 Feature f -> (Settings (gtf settings)
                                                                        (fasta settings)
                                                                        f
                                                                        (length' settings)
                                                                        (offset settings)
                                                                        (phase' settings), count )
                                                 Phase p -> (Settings (gtf settings)
                                                                      (fasta settings)
                                                                      (feature settings)
                                                                      (length' settings)
                                                                      (offset settings)
                                                                      ((read (n!!count)) :: Integer), count + 1)

             defaults = Settings "" "" ""  9 3 (-1)




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
                   | feature s == "all-exons" = extractExons s a
                   | feature s == "intergenic" = extractIntergenic s a
                   | feature s == "single" = extractSingleExons s a
                   | feature s == "cds" = extractCDS s a
                   | otherwise = extractExons s a



extractIntergenic :: Settings -> [Annotation] -> ([String], [String])
extractIntergenic s a = (concat $ map (\ annot -> printIntergenic annot) a, [])
                  where printIntergenic annot = map (\ positions -> printSequence name (fst positions) (snd positions) ) (getIntergenicRegions annot)
                          where name = seqname (seqentry annot)



extractIntrons :: Settings -> [Annotation] -> ([String], [String])
extractIntrons s a = extractFeature s a fromGenes (\ annot -> genes annot)
                    where fromGenes s p g = extractFeature s g fromTranscripts ( \ gene -> transcripts gene )
                          fromTranscripts s p t = extractFeature s t fromIntrons ( \ txs -> (getIntrons txs) )
                          fromIntrons s p introns = case (strand p == "+") of
                                                       True ->  (seqstr, [])
                                                       False -> ([],seqstr)
                                                    where seqstr = map ( \ i -> printSequence seqname (fst i) (snd i) ) introns
                                                          seqname = parentseq (cstart ((txcds p) !! 0))



extractSingleExons :: Settings -> [Annotation] -> ([String], [String])
extractSingleExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts ( \ x -> transcripts x )
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
                                                   seqstr = printSequence (parentseq (cstart c))
                                                                          (position (cstart c))
                                                                          (position (cend c))



extractFinalExons :: Settings -> [Annotation] -> ([String], [String])
extractFinalExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts ( \ x -> transcripts x )
                       fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithStartPhase (phase' s) x) )
                       fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                       fromCDS s p (c:_)= case (stype (cstart c)) of
                                             Acceptor -> case (stype (cend c)) of
                                                           StopCodon -> forward
                                                           _ -> ([],[])
                                             StopCodon -> case (stype (cend c)) of
                                                           Acceptor -> reverse
                                                           _ -> ([],[])
                                             _ -> ([],[])
                                             where forward = ([seqstr],[])
                                                   reverse = ([], [seqstr])
                                                   seqstr = printSequence (parentseq (cstart c))
                                                                          (position (cstart c))
                                                                          (position (cend c))



extractInternalExons :: Settings -> [Annotation] -> ([String], [String])
extractInternalExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts ( \ x -> transcripts x )
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
                                                   seqstr = printSequence (parentseq (cstart c))
                                                                          (position (cstart c))
                                                                          (position (cend c))



extractInitialExons :: Settings -> [Annotation] -> ([String], [String])
extractInitialExons s a = extractFeature s a fromGenes (\ x -> genes x)
                 where fromGenes s p g = extractFeature s g fromTranscripts ( \ x -> transcripts x )
                       fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithEndPhase (phase' s) x))
                       fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                       fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> forward
                                             Donor -> case (stype (cend c)) of
                                                        StartCodon -> reverse
                                                        _ -> ([],[])
                                             _ -> ([],[])
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstr])
                                                  seqstr = printSequence (parentseq (cstart c))
                                                                         (position (cstart c))
                                                                         (position (cend c))


extractCDS :: Settings -> [Annotation] -> ([String], [String])
extractCDS s a = extractFeature s a fromGenes (\ x -> genes x)
                   where fromGenes s p g = extractFeature s g fromTranscripts ( \ x -> transcripts x )
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
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstr])
                                                  seqstr = " "++ show (position (cstart c)) ++ " " ++ show  (position (cend c))



extractExons :: Settings -> [Annotation] -> ([String], [String])
extractExons s a = extractFeature s a fromGenes (\ x -> genes x)
                   where fromGenes s p g = extractFeature s g fromTranscripts ( \ x -> transcripts x )
                         fromTranscripts s p t = extractFeature s t fromCDSList ( \x -> (getCDSWithStartPhase (phase' s) x))
                         fromCDSList s p c = extractFeature s c fromCDS ( \x -> [x] )
                         fromCDS s p (c:_)= case (stype (cstart c)) of
                                             StartCodon -> forward
                                             Acceptor -> forward
                                             Donor -> reverse
                                             StopCodon -> reverse
                                            where forward = ([seqstr],[])
                                                  reverse = ([], [seqstr])
                                                  seqstr = printSequence (parentseq (cstart c))
                                                                         (position (cstart c))
                                                                         (position (cend c))


extractSites :: Settings -> [Annotation] -> ([String], [String])
extractSites s a = extractFeature s a fromGenes (\ annot -> genes annot)
                    where fromGenes s p g = extractFeature s g fromTranscripts ( \ gene -> transcripts gene )
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
printSequence n s e = n ++ " " ++ " " ++ show (s) ++ " " ++ show e ++ "\n"

