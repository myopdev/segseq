
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
data Transcript = Transcript { txname::Name, strand::Strand, startSite::Site, txcds::[CDS] , endSite::Site }
                  deriving (Show)
data Gene = Gene {gname :: Name , transcripts::[Transcript]}
            deriving (Show)
data Annotation = Annotation {seqentry :: Sequence, genes :: [Gene], source :: String}
                  deriving (Show)


sequence2 = Sequence "TEST2" 0 1000
site5 =  Site "TEST2" 1 StartCodon
site6 =  Site "TEST2" 100 StopCodon
site7 = Site "TEST2" 1 StopCodon
site8 = Site "TEST2" 100 StartCodon
cds2 = CDS site5 site6 0
tx4 = Transcript "C.1" "+" site5 [cds2] site6
tx5 = Transcript "D.1" "-" site7 [cds2] site8
tx6 = Transcript "D.2" "+" site5 [cds2] site6
gene3 = Gene "C" [tx4]
gene4 = Gene "D" [tx5, tx6]
annot1 = Annotation sequence2 [gene3, gene4] "tops"


sequence1 = Sequence "TEST1" 0 1000
site1 =  Site "TEST1" 1 StartCodon
site2 =  Site "TEST1" 100 StopCodon
site3 = Site "TEST1" 1 StopCodon
site4 = Site "TEST1" 100 StartCodon
cds1 = CDS site1 site2 0
tx1 = Transcript "A.1" "+" site1 [cds1] site2
tx2 = Transcript "B.1" "-" site3 [cds1] site4
tx3 = Transcript "B.2" "+" site1 [cds1] site2
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




