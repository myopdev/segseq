module SegSeq.SequenceAnnotation where

import Data.List
import Maybe

type Name = String
type Strand = String

data SiteType = StartCodon
              | StopCodon
              | Acceptor
              | Donor
              | Unknown
                deriving (Show)
data Sequence = Sequence {seqname::Name, seqstart:: Integer, seqend::Integer}
                deriving (Show)
data Site = Site {parentseq::Name, position:: Integer, stype:: SiteType}
            deriving (Show)
data CDS = CDS {cstart::Site, cend::Site, phase::Integer}
           deriving (Show)
data Transcript = Transcript { txname::Name, strand::Strand, txcds::[CDS], gcAll::Integer, gcFromCDS::Integer, gcFromIntron::Integer }
                  deriving (Show)
data Gene = Gene {gname :: Name , transcripts::[Transcript]}
            deriving (Show)
data Annotation = Annotation {seqentry :: Sequence, genes :: [Gene], source :: String}
                  deriving (Show)


getgeneids :: Annotation -> [String]
getgeneids a = Prelude.map ( \ x -> gname x) (genes a)

gettranscriptids :: Gene -> [String]
gettranscriptids g = Prelude.map ( \x -> txname x ) (transcripts g)

gettranscript :: Annotation -> String -> Transcript
gettranscript a id = fromJust (find ( \ x -> ( (txname x) == id ) ) txs)
              where txs = concat (Prelude.map ( \x -> transcripts x ) (genes a))


getCDSWithStartPhase ::   Integer -> Transcript -> [CDS]
getCDSWithStartPhase phase' tx  | phase' >= 0  =  filter (\ x -> (phase x == phase')) (txcds tx)
                                | otherwise = txcds tx

getCDSWithEndPhase ::  Integer -> Transcript -> [CDS]
getCDSWithEndPhase phase' tx | phase' >= 0  =  filter (\ x -> (((phase x) + (l x) - 1) `mod` 3) == (phase' )  ) (txcds tx)
                             | otherwise = txcds tx
                            where l x= ((position $cend x) - (position $cstart x)+1)

getgene :: Annotation -> String -> Gene
getgene a id = fromJust (find ( \ x ->( (gname x) == id) ) (genes a))


getSites :: Integer-> SiteType ->Transcript -> [Site]
getSites phase' stype' = case stype' of
                        StartCodon -> getStartCodonSites phase'
                        StopCodon -> getStopCodonSites phase'
                        Acceptor -> getAcceptorSites phase'
                        Donor -> getDonorSites phase'

getAcceptorSites :: Integer -> Transcript -> [Site]
getAcceptorSites phase' tx  =
                     case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cstart cds)  (getCDSWithStartPhase  phase' tx ))
                       False -> filter isSite (map (\ cds -> cend cds)  (getCDSWithStartPhase phase' tx))
                     where isSite site =  case stype site of
                                               Acceptor -> True
                                               _ -> False

getStartCodonSites :: Integer -> Transcript  -> [Site]
getStartCodonSites phase' tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cstart cds)  (getCDSWithStartPhase  phase' tx))
                       False -> filter isSite (map (\ cds -> cend cds)  (getCDSWithStartPhase  phase' tx))
                     where isSite site =  case stype site of
                                               StartCodon -> True
                                               _ -> False

getDonorSites :: Integer -> Transcript  -> [Site]
getDonorSites phase' tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cend cds)  (getCDSWithEndPhase  phase' tx))
                       False -> filter isSite (map (\ cds -> cstart cds)  (getCDSWithEndPhase  phase' tx))
                     where isSite site =  case stype site of
                                               Donor -> True
                                               _ -> False

getStopCodonSites :: Integer -> Transcript  -> [Site]
getStopCodonSites phase' tx = case (strand tx) == "+" of
                       True -> filter isSite (map (\ cds -> cend cds)  (getCDSWithEndPhase  phase' tx))
                       False -> filter isSite (map (\ cds -> cstart cds)  (getCDSWithEndPhase  phase' tx))
                     where isSite site =  case stype site of
                                               StopCodon -> True
                                               _ -> False

getIntrons :: Transcript -> [(Integer, Integer)]
getIntrons tx = case (strand tx ) == "+" of
                 True -> getIntronForward tx
                 False -> getIntronReverse tx
               where getIntronForward tx = zip  (map (addPosition 1)  $ getDonorSites (-1) tx) (map (addPosition (-1)) $ getAcceptorSites (-1) tx)
                     getIntronReverse tx = zip  (map (addPosition 1)  $ getAcceptorSites (-1) tx) (map (addPosition (-1)) $ getDonorSites (-1) tx)
                     addPosition n site  = ((position site) +n)

getExons :: Integer->Transcript -> [(Integer, Integer)]
getExons phase' tx = map (\ cds -> ((position $ cstart cds), (position $ cend cds))) (getCDSWithStartPhase phase' tx)

getGenePosition :: Gene -> (Integer, Integer)
getGenePosition g = (minimum positions, maximum positions)
           where positions = concat $ map (\ cds -> [(position $ cstart cds), (position $cend cds)] ) ( concat $ map (\ tx -> txcds tx) (transcripts g))

getTranscriptPosition :: Transcript -> (Integer, Integer)
getTranscriptPosition tx = (minimum positions, maximum positions)
  where positions = concat $ map (\ cds -> [(position $ cstart cds), (position $cend cds)] ) ( txcds tx)

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
                                       complement [(a,b)] = []
                                       complement ([(a,b),(c,d)]) = case (b < c - 1) of
                                                                      True -> [(b+1, c-1), (d+1, d+1000)]
                                                                      False ->[(d+1, d+1000)]
                                       complement ((a,b):(c,d):xs) =  [(b+1, c-1)] ++ (complement ((c,d):xs))


getIntergenicRegions :: Annotation -> [(Integer, Integer)]
getIntergenicRegions a = complementIntervals $ disjointIntervals  (map (\ gene -> getGenePosition gene ) (genes a))


getNonCodingRegions :: Annotation -> [(Integer, Integer)]
getNonCodingRegions a = complementIntervals $ disjointIntervals  $ concat $ map (\ gene -> (concat $ map (getExons (-1))  (transcripts  gene) )) (genes a)



