module SequenceAnnotation where

import Data.List
import Maybe

type Name = String
type Strand = String

data SiteType = StartCodon
              | StopCodon
              | Acceptor
              | Donor
                deriving (Show)
data Sequence = Sequence {seqname::Name, seqstart:: Integer, seqend::Integer}
                deriving (Show)
data Site = Site {parentseq::Name, position:: Integer, stype:: SiteType}
            deriving (Show)
data CDS = CDS {cstart::Site, cend::Site, phase::Integer}
           deriving (Show)
data Transcript = Transcript { txname::Name, strand::Strand, txcds::[CDS] }
                  deriving (Show)
data Gene = Gene {gname :: Name , transcripts::[Transcript]}
            deriving (Show)
data Annotation = Annotation {seqentry :: Sequence, genes :: [Gene], source :: String}
                  deriving (Show)


sequence2 = Sequence "TEST2" 0 1000
site5 =  Site "TEST2" 1 StartCodon
site6 =  Site "TEST2" 100 Donor
site9 =  Site "TEST2" 120 Acceptor
site10 =  Site "TEST2" 130 Donor
site11 =  Site "TEST2" 150 Acceptor
site12 =  Site "TEST2" 160 StopCodon
site7 = Site "TEST2" 1 StopCodon
site8 = Site "TEST2" 100 StartCodon
cds2 = CDS site5 site6 0
cds3 = CDS site9 site10 0
cds4 = CDS site11 site12 0
tx4 = Transcript "C.1" "+"  [cds2, cds3, cds4]
gene3 = Gene "C" [tx4]
annot1 = Annotation sequence2 [gene3] "tops"


sequence1 = Sequence "TEST1" 0 1000
site1 =  Site "TEST1" 100 StartCodon
site2 =  Site "TEST1" 140 StopCodon
site3 = Site "TEST1" 110 StopCodon
site4 = Site "TEST1" 120 StartCodon
cds1 = CDS site1 site2 0
cds5 = CDS site3 site4 0
tx1 = Transcript "A.1" "+"  [cds1]
tx2 = Transcript "B.1" "-"  [cds5]
tx3 = Transcript "B.2" "+"  [cds1]
gene1 = Gene "A" [tx1]
gene2 = Gene "B" [tx2, tx3]
annot2 = Annotation sequence1 [gene1, gene2] "test"



getgeneids :: Annotation -> [String]
getgeneids a = Prelude.map ( \ x -> gname x) (genes a)

gettranscriptids :: Gene -> [String]
gettranscriptids g = Prelude.map ( \x -> txname x ) (transcripts g)

gettranscript :: Annotation -> String -> Transcript
gettranscript a id = fromJust (find ( \ x -> ( (txname x) == id ) ) txs)
              where txs = concat (Prelude.map ( \x -> transcripts x ) (genes a))


getcds :: Transcript -> [CDS]
getcds tx = txcds tx

getgene :: Annotation -> String -> Gene
getgene a id = fromJust (find ( \ x ->( (gname x) == id) ) (genes a))


getSites :: SiteType -> Transcript -> [Site]
getSites stype' = case stype' of
                        StartCodon -> getStartCodonSites
                        StopCodon -> getStopCodonSites
                        Acceptor -> getAcceptorSites
                        Donor -> getDonorSites

getAcceptorSites :: Transcript  -> [Site]
getAcceptorSites tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cstart cds)  (txcds tx))
                       False -> filter isSite (map (\ cds -> cend cds)  (txcds tx))
                     where isSite site =  case stype site of
                                               Acceptor -> True
                                               _ -> False

getStartCodonSites :: Transcript  -> [Site]
getStartCodonSites tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cstart cds)  (txcds tx))
                       False -> filter isSite (map (\ cds -> cend cds)  (txcds tx))
                     where isSite site =  case stype site of
                                               StartCodon -> True
                                               _ -> False

getDonorSites :: Transcript  -> [Site]
getDonorSites tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cend cds)  (txcds tx))
                       False -> filter isSite (map (\ cds -> cstart cds)  (txcds tx))
                     where isSite site =  case stype site of
                                               Donor -> True
                                               _ -> False

getStopCodonSites :: Transcript  -> [Site]
getStopCodonSites tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cend cds)  (txcds tx))
                       False -> filter isSite (map (\ cds -> cstart cds)  (txcds tx))
                     where isSite site =  case stype site of
                                               StopCodon -> True
                                               _ -> False

getIntrons :: Transcript -> [(Integer, Integer)]
getIntrons tx = case (strand tx ) == "+" of
                 True -> getIntronForward tx
                 False -> getIntronReverse tx
               where getIntronForward tx = zip  (map (addPosition 1)  $ getDonorSites tx) (map (addPosition (-1)) $ getAcceptorSites tx)
                     getIntronReverse tx = zip  (map (addPosition 1)  $ getAcceptorSites tx) (map (addPosition (-1)) $ getDonorSites tx)
                     addPosition n site  = ((position site) +n)

getExons :: Transcript -> [(Integer, Integer)]
getExons tx = map (\ cds -> ((position $ cstart cds), (position $ cend cds))) (txcds tx)

getGenePosition :: Gene -> (Integer, Integer)
getGenePosition g = (minimum positions, maximum positions)
           where positions = concat $ map (\ cds -> [(position $ cstart cds), (position $cend cds)] ) ( concat $ map (\ tx -> txcds tx) (transcripts g))

disjointIntervals :: [(Integer, Integer)] -> [(Integer, Integer)]
disjointIntervals intervals = joinIntervals (sort intervals)
                  where joinIntervals [] = []
                        joinIntervals ([(a,b)]) = [(a,b)]
                        joinIntervals ([(a,b),(c,d)]) =  case ((a <= c ) && (c <= b+1)) of
                                                           True -> [(a,e)] where e = maximum [b,d]
                                                           False -> [(a,b), (c,d)]
                        joinIntervals ((a,b):(c,d):xs) = case ((a <= c ) && (c <= b+1)) of
                                                           True -> joinIntervals ((a,e):xs) where e = maximum [b,d]
                                                           False -> [(a,b)] ++ joinIntervals ((c,d):xs)

complementIntervals :: [(Integer, Integer)] -> [(Integer, Integer)]
complementIntervals intervals = complement ((0,0):intervals)
                                 where complement [] = []
                                       complement ([(a,b),(c,d)]) = case (b < c - 1) of
                                                                      True -> [(b+1, c-1), (d+1, d+1000)]
                                                                      False ->[(d+1, d+1000)]
                                       complement ((a,b):(c,d):xs) =  [(b+1, c-1)] ++ (complement ((c,d):xs))


getIntergenicRegions :: Annotation -> [(Integer, Integer)]
getIntergenicRegions a = complementIntervals $ disjointIntervals  (map (\ gene -> getGenePosition gene ) (genes a))


