module GTF where
import SequenceAnnotation
import Text.ParserCombinators.Parsec hiding (spaces)
import Data.CSV

renderGTF :: Annotation -> String
renderGTF  a = genesGTF a (genes a)

genesGTF :: Annotation -> [Gene] -> String
genesGTF a [] = ""
genesGTF a (x:xs) = geneGTF a x  ++ genesGTF a xs

geneGTF :: Annotation -> Gene -> String
geneGTF a g = transcriptsGTF a g (transcripts g)

transcriptsGTF :: Annotation -> Gene -> [Transcript] -> String
transcriptsGTF a g [] = ""
transcriptsGTF a g (x:xs) = transcriptGTF a g x ++  transcriptsGTF a g xs

transcriptGTF :: Annotation -> Gene -> Transcript -> String
transcriptGTF a g t = cdsListGTF a g t (txcds t)

cdsListGTF :: Annotation -> Gene -> Transcript -> [CDS] -> String
cdsListGTF a g t [] = ""
cdsListGTF a g t (x:xs) = (siteGTF a g t (startSite t)) ++ "\n" ++  (cdsGTF a g t x) ++  "\n" ++ (siteGTF a g t (endSite t)) ++ "\n" ++ (cdsListGTF a g t xs)

siteGTF :: Annotation -> Gene-> Transcript -> Site -> String
siteGTF a g t s = case (strand t ) of "+" -> forwardStrandSiteGTF a g t s
                                      "-" -> reverseStrandSiteGTF a g t s

forwardStrandSiteGTF :: Annotation -> Gene-> Transcript -> Site -> String
forwardStrandSiteGTF a g t s = case (stype s) of StartCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "start_codon" ++ "\t" ++ (show (position s)) ++ "\t" ++ (show ((position s) + 2)) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\""
                                                 StopCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "stop_codon" ++ "\t" ++ (show ((position s) + 1)) ++ "\t" ++ (show ((position s) + 3)) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\"\n"
                                                 _ -> ""

reverseStrandSiteGTF :: Annotation -> Gene-> Transcript -> Site -> String
reverseStrandSiteGTF a g t s =  case (stype s) of StartCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "start_codon" ++ "\t" ++ (show ((position s) - 2)) ++ "\t" ++ (show ((position s) )) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\"\n"
                                                  StopCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "stop_codon" ++ "\t" ++ (show ((position s) - 3)) ++ "\t" ++ (show ((position s) - 1)) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\""
                                                  _ -> ""

cdsGTF ::  Annotation -> Gene -> Transcript -> CDS -> String
cdsGTF a g t c = seqname1 ++ "\t" ++ source1 ++ "\t" ++ "CDS" ++ "\t" ++ start1 ++ "\t" ++ stop1 ++ "\t" ++ ".\t" ++  strand1 ++ "\t" ++ phase1 ++ "\tgene_id \"" ++ gname1 ++ "\"; "++ "transcript_id \"" ++ txname1 ++ "\""
              where
                    seqname1 = seqname (seqentry a)
                    source1 = source a
                    start1 = show (position (cstart c))
                    stop1 =  show (position (cend c))
                    phase1 = show (phase c)
                    strand1 = strand t
                    gname1 = gname g
                    txname1 = txname t


gtfparser :: GenParser Char st [[String]]
gtfparser = do result <- many line
               eof
               return result
line :: GenParser Char st [String]
line = do result <- cells
          eol
          return result
cells:: GenParser Char st [String]
cells =  do first <- cellContent
            next <- remainingCells
            return (first : next)

remainingCells :: GenParser Char st [String]
remainingCells = (char '\t' >> cells ) <|> (return [])

cellContent :: GenParser Char st String
cellContent = many (noneOf "\t\n")

eol :: GenParser Char st Char
eol = char '\n'

readCSV input = case parse gtfparser "GTF" input of
        Left err -> [["Error: " ++ show err]]
        Right val -> val


-- readGTF :: String -> Annotation
-- readGTF x = (initializeAnnotationGTF . readCSV) x






