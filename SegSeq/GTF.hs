module SegSeq.GTF where
import SegSeq.SequenceAnnotation
import Text.ParserCombinators.Parsec
import Data.List

renderGTF :: Annotation -> String
renderGTF a = renderGenesGTF a (genes a)

renderGenesGTF :: Annotation -> [Gene] -> String
renderGenesGTF a [] = ""
renderGenesGTF a (x:xs) = renderGeneGTF a x  ++ renderGenesGTF a xs

renderGeneGTF :: Annotation -> Gene -> String
renderGeneGTF a g = renderTranscriptsGTF a g (transcripts g)

renderTranscriptsGTF :: Annotation -> Gene -> [Transcript] -> String
renderTranscriptsGTF a g [] = ""
renderTranscriptsGTF a g (x:xs) = renderTranscriptGTF a g x ++  renderTranscriptsGTF a g xs

renderTranscriptGTF :: Annotation -> Gene -> Transcript -> String
renderTranscriptGTF a g t = (renderCDSListGTF a g t (txcds t))

renderCDSListGTF :: Annotation -> Gene -> Transcript -> [CDS] -> String
renderCDSListGTF a g t [] = ""
renderCDSListGTF a g t (x:xs) =    (renderCDSGTF a g t x)  ++ (renderCDSListGTF a g t xs)


renderCDSGTF ::  Annotation -> Gene -> Transcript -> CDS -> String
renderCDSGTF a g t c =  begin ++  seqname1 ++ "\t" ++ source1 ++ "\t" ++ "CDS" ++ "\t" ++ start1 ++ "\t" ++ stop1 ++ "\t" ++ ".\t" ++  strand1 ++ "\t" ++ phase1 ++ "\tgene_id \"" ++ gname1 ++ "\"; "++ "transcript_id \"" ++ txname1 ++ "\";\n" ++  end
                        where begin =  renderSiteGTF a g t type1
                              end = renderSiteGTF a g t type2
                              seqname1 = seqname (seqentry a)
                              source1 = source a
                              start1 = show (position (cstart c))
                              stop1 =  show (position (cend c))
                              phase1 = show (phase c)
                              strand1 = strand t
                              gname1 = gname g
                              txname1 = txname t
                              type1 =   cstart c
                              type2 =  cend c

renderSiteGTF :: Annotation -> Gene-> Transcript -> Site -> String
renderSiteGTF a g t s = case (strand t ) of "+" -> forwardStrandSiteGTF a g t s
                                            "-" -> reverseStrandSiteGTF a g t s


forwardStrandSiteGTF :: Annotation -> Gene-> Transcript -> Site -> String
forwardStrandSiteGTF a g t s = case (stype s) of StartCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "start_codon" ++ "\t" ++ (show (position s)) ++ "\t" ++ (show ((position s) + 2)) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\";\n"
                                                 StopCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "stop_codon" ++ "\t" ++ (show ((position s) + 1)) ++ "\t" ++ (show ((position s) + 3)) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\";\n"
                                                 _ -> ""

reverseStrandSiteGTF :: Annotation -> Gene-> Transcript -> Site -> String
reverseStrandSiteGTF a g t s =  case (stype s) of StartCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "start_codon" ++ "\t" ++ (show ((position s) - 2)) ++ "\t" ++ (show ((position s) )) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\";\n"
                                                  StopCodon -> (seqname (seqentry a)) ++ "\t" ++ (source a)++ "\t" ++ "stop_codon" ++ "\t" ++ (show ((position s) - 3)) ++ "\t" ++ (show ((position s) - 1)) ++ "\t" ++ ".\t" ++  strand t ++ "\t0\tgene_id \"" ++ (gname g) ++ "\"; "++ "transcript_id \"" ++ (txname t) ++ "\";\n"
                                                  _ -> ""



gtfparser :: GenParser Char st  [[String]]
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

readCSV :: String ->  [[String]]
readCSV input = case parse gtfparser "GTF" input of
                     Left err -> [["Parse error: " ++ show err]]
                     Right val ->  removeEmptyEntries val

se :: Bool -> a -> a -> a
se True  x _ = x
se False _ y = y

removeEmptyEntries ::[[String]] ->[[String]]
removeEmptyEntries []  = []
removeEmptyEntries (x:xs) = se (x == [""]) (removeEmptyEntries xs) ([x] ++ removeEmptyEntries xs)

readAttr :: String ->  [String]
readAttr input = case parse (attrParser <|> attrParserR) "parsing attribute" input of
                      Left err -> ["Parse error: " ++ show err]
                      Right val -> val

p_string :: GenParser Char st String
p_string = do char '"'
              t <- many (noneOf "\"")
              char '"'
              return t

attrParser :: GenParser Char st [String]
attrParser  = do geneid <- geneIdParser
                 transcriptid <- txParser
                 return [geneid, transcriptid];


attrParserR :: GenParser Char st [String]
attrParserR  = do transcriptid <- txParser
                  geneid <- geneIdParser
                  return [geneid, transcriptid];


geneIdParser :: GenParser Char st String
geneIdParser = do  string "gene_id" >> spaces
                   g <- p_string
                   string ";" >> spaces
                   return g


txParser :: GenParser Char st String
txParser = do string "transcript_id" >> spaces
              t <- p_string
              string ";" >> spaces
              return t


buildGeneNameList :: [[String]] -> [String]
buildGeneNameList gtflines = nub $ map getGeneNameFromGTFLine  (delete [""] (nub gtflines))

buildTranscriptNameList :: [[String]] -> [String]
buildTranscriptNameList gtflines = nub $ map getTranscriptNameFromGTFLine  (delete [""] (nub gtflines))


clusterGTFLinesBySequence :: [[String]] -> [[[String]]]
clusterGTFLinesBySequence gtflines = map (\ x -> filter ( \ y -> ((y !! 0) ==  (x!!0)) && ((y !! 1) ==  (x!!1))  ) gtflines )  $ sequenceNames
                                     where sequenceNames = nub $ map (\ x -> ([x !! 0, x!!1])) $ delete [""] (nub gtflines)


createAnnotationList :: [[String]] -> [Annotation]
createAnnotationList gtflines = map createAnnotation $ clusterGTFLinesBySequence gtflines

createAnnotation ::[[String]] -> Annotation
createAnnotation gtflines = Annotation seq genes source
                                       where seq = Sequence ((gtflines !! 0) !! 0) 0 (maximum (map (\x -> read $x !! 3) gtflines ))
                                             genes = createGeneList gtflines
                                             source = gtflines !! 0  !! 1

sortByStartPosition :: [[String]] -> [[String]]
sortByStartPosition gtflines = sortBy comparaGTFLine gtflines

comparaGTFLine :: [String] -> [String] -> Ordering
comparaGTFLine x y
  | ((read ( x!!3)::Integer ) > (read (y!!3 )::Integer)) = GT
  | otherwise = LT




createGeneList :: [[String]] -> [Gene]
createGeneList gtflines = map createGene $ clusterGTFLinesByGeneName gtflines


createGene :: [[String]] -> Gene
createGene gtflines = Gene name transcripts
                           where name = (getGeneNameFromGTFLine (gtflines !! 0))
                                 transcripts = map createTranscript $ clusterGTFLinesByTranscript gtflines

createTranscript :: [[String]] -> Transcript
createTranscript gtflines = Transcript name strand cds 0 0 0
                                       where name = getTranscriptNameFromGTFLine (gtflines !! 0)
                                             strand = (gtflines !! 0) !! 6
                                             cds = calculateCDSPhase (sortBy compareCDS ((createCDSList gtflines (createStartSite gtflines))  (createStopSite gtflines))) strand

compareCDS :: CDS -> CDS -> Ordering
compareCDS x y | (position $ cstart x) > (position $ cstart y)  = GT
               | otherwise = LT

calculateCDSPhase :: [CDS] -> String -> [CDS]
calculateCDSPhase lcds strand = se (strand == "+")
                                  (foldl (\ listCDS  c -> listCDS ++ [CDS (cstart c) (cend c) (getPhaseForward listCDS) ])  [head lcds] (tail  lcds))
                                  (foldr (\ c listCDS -> [CDS (cstart c) (cend c) (getPhaseReverse listCDS) ] ++ listCDS)  [last lcds] (init lcds))
                             where getPhaseForward l = (((position $cend $last(l)) - (position $cstart $last(l)) + 1) + phase (last (l))) `mod` 3
                                   getPhaseReverse l = (((position $cend $head(l)) - (position $cstart $head(l)) + 1) + phase (head (l))) `mod` 3

createStartSite :: [[String]] -> Site
createStartSite gtflines = se (length y > 0)
                              (se ((y !! 0 !! 6) == "+" )
                                  (Site (y!! 0 !! 0) (read (y !! 0 !! 3)) StartCodon)
                                  (Site (y!! 0 !! 0) (read (y !! 0 !! 4)) StartCodon) )
                              (Site "Error" (-1000) Unknown)
                           where  y = ( filter ( \ x -> ( (x !! 2) == "start_codon" )) gtflines )

createStopSite :: [[String]] -> Site
createStopSite gtflines =  se (length y > 0)
                              (se ((y !! 0 !! 6) == "+" )
                                  (Site (y!! 0 !! 0) (read (y !! 0 !! 3) -1) StopCodon)
                                  (Site (y!! 0 !! 0) (read(y !! 0 !! 4)+1) StopCodon))
                              (Site "Error" (-1000) Unknown)
                           where  y = ( filter ( \ x -> ( (x !! 2) == "stop_codon" )) gtflines )


createCDSList :: [[String]] -> Site-> Site -> [CDS]
createCDSList gtflines site1 site2 = map ( \ x -> createCDS x site1 site2)  ( filter ( \ x -> ( (x !! 2) == "CDS" )) gtflines )

createCDS :: [String] -> Site -> Site -> CDS
createCDS gtfline site1 site2
          | (hasUnknown site1 site2) = CDS (Site name pos1 (stype site1)) (Site name pos2 (stype site2)) phase
          | (position site1 == pos1 && position site2 == pos2) = CDS (Site name pos1 StartCodon) (Site name pos2 StopCodon) phase
          | (position site1 == pos2 && position site2 == pos1) = CDS (Site name pos1 StopCodon) (Site name pos2 StartCodon) phase
          | (position site1 == pos1 && position site2 /= pos2) = CDS (Site name pos1 StartCodon) (Site name pos2 Donor) phase
          | (position site1 == pos2 && position site2 /= pos1) = CDS (Site name pos1 Donor) (Site name pos2 StartCodon) phase
          | (position site1 /= pos1 && position site2 == pos2) = CDS (Site name pos1 Acceptor) (Site name pos2 StopCodon) phase
          | (position site1 /= pos2 && position site2 == pos1) = CDS (Site name pos1 StopCodon) (Site name pos2 Acceptor) phase
          | otherwise = se (strand == "+") (CDS (Site name pos1 Acceptor) (Site name pos2 Donor) phase)
                                           (CDS (Site name pos1 Donor) (Site name pos2 Acceptor) phase)
          where pos1 = read $ gtfline !! 3
                pos2 = read $ gtfline !! 4
                strand = gtfline !! 6
                name = gtfline !! 0
                phase = read $ gtfline !! 7
                hasUnknown s1 s2 = case (stype site1) of
                                      Unknown -> True
                                      _ -> case (stype site2) of
                                              Unknown -> True
                                              _ -> False

clusterGTFLinesByTranscript :: [[String]] ->[[[String]]]
clusterGTFLinesByTranscript gtflines = map (\ x -> getTranscript x gtflines ) $ buildTranscriptNameList gtflines

clusterGTFLinesByGeneName :: [[String]] -> [[[String]]]
clusterGTFLinesByGeneName gtflines = map (\ x -> getGene x gtflines)  $ buildGeneNameList gtflines

getGene :: String -> [[String]] -> [[String]]
getGene name gtflines = filter (\x -> (getGeneNameFromGTFLine x) == name ) gtflines

getTranscript :: String -> [[String]] -> [[String]]
getTranscript name gtflines = filter (\x -> (getTranscriptNameFromGTFLine x) == name ) gtflines

getGeneNameFromGTFLine :: [String] -> String
getGeneNameFromGTFLine s = case length s of
                                9 -> (readAttr $ s !! 8) !! 0
                                _ -> ""

getTranscriptNameFromGTFLine :: [String] -> String
getTranscriptNameFromGTFLine s = case length s of
                                9 -> (readAttr $ s !! 8) !! 1
                                _ -> ""

