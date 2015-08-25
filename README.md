# segseq

A simple program that receives a GTF (Gene Transfer Format) file and a FASTA file and returns the sequences that matches with 
a specific annotated feature.

# Software requirement

1. [Haskell](https://www.haskell.org/)
2. [Cabal](https://www.haskell.org/cabal)

# Installing SegSeq

1. Download *segseq* from GitHub 

```
git clone https://github.com/myopdev/segseq.git
```

2. Execute *cabal*

```
cabal install segseq
```


# Examples 

## Acceptor splice sites

``` 
segseq -r -x acceptor -f dataset/genes.fa -g dataset/genes.gtf -l 6 -o 3 > acceptor.fa
```

## Donor splice sites 
``` 
segseq -r -x donor -l 9 -o 3 -f dataset/genes.fa -g dataset/genes.gtf > donor.fa
```


