
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
site1 =  Site "TEST1" 1 StartCodon
site2 =  Site "TEST1" 100 StopCodon
site3 = Site "TEST1" 1 StopCodon
site4 = Site "TEST1" 100 StartCodon
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




