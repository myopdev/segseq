
module SequenceAnnotation where

import Data.List


type Name = String
type Strand = String
data SiteType = StartCodon
              | StopCodon
              | Acceptor
              | Donor
                deriving (Show)

data Sequence = Sequence {seqname::Name, seqstart:: Integer, seqend::Integer}
                deriving (Show)
data Site = Site {position:: Integer, stype:: SiteType}
            deriving (Show)
data CDS = CDS {cstart::Site, cend::Site, phase::Integer}
           deriving (Show)
data Transcript = Transcript { txname::Name, strand::Strand, startSite::Site, txcds::[CDS] , endSite::Site }
                  deriving (Show)
data Gene = Gene {gname :: Name , transcripts::[Transcript]}
            deriving (Show)
data Annotation = Annotation {seqentry :: Sequence, genes :: [Gene], source :: String}
                  deriving (Show)

data AnnotationVal = AnnotationVal Annotation
                   | GeneVal Gene
                   | TranscriptVal Transcript
                   | CDSVal CDS
                   | SiteVal Site
                   | SequenceVal Sequence

sequence1 = Sequence "TEST1" 0 1000
site1 = Site 1 StartCodon
site2 = Site 100 StopCodon
site3 = Site 1 StopCodon
site4 = Site 100 StartCodon

cds1 = CDS site1 site2 0
tx1 = Transcript "A.1" "+" site1 [cds1] site2
tx2 = Transcript "B.1" "-" site3 [cds1] site4
tx3 = Transcript "B.2" "+" site1 [cds1] site2
gene1 = Gene "A" [tx1]
gene2 = Gene "B" [tx2, tx3]
annot = Annotation sequence1 [gene1, gene2] "test"


geneids :: Annotation -> [String]
geneids a = Prelude.map ( \ x -> gname x) (genes a)

transcriptids :: Gene -> [String]
transcriptids g = Prelude.map ( \x -> txname x ) (transcripts g)

getgene :: Annotation -> String -> Maybe Gene
getgene a id = find ( \ x ->( (gname x) == id) ) (genes a)



