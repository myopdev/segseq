# segseq

A simple program that receives a GTF (Gene Transfer Format) file and a FASTA file and returns the sequences that matches with 
a specific annotated feature.

# Software requirement

1. [Haskell](https://www.haskell.org/)
2. [Cabal](https://www.haskell.org/cabal)
3. [SQLite3](https://www.sqlite.org/):  ``libsqlite3-dev``

# Installing SegSeq

1. Download *segseq* from GitHub 

 ```
 git clone https://github.com/myopdev/segseq.git
 ```

2. Execute *cabal update*

 ``` 
 cabal update
 ```

3. Execute *cabal install*

 ```
 cabal install ./segseq
 ```
4. Add the following line at the end of the ``.profile`` file

 ```
 export PATH=$PATH:$HOME/.cabal/bin
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

### CDS

```
segseq -r -x cds -f dataset/genes.fa -g dataset/genes.gtf
```

### Intron

```
 segseq  -r -x intron -g dataset/genes.gtf -f dataset/genes.fa >
```

### Initial Exons

```
segseq -r -x initial -g dataset/genes.gtf -f dataset/genes.fa
```

### Internal Exons

```
segseq -r -x internal -g dataset/genes.gtf -f dataset/genes.fa
```

### Final Exons

```
segseq -r -x final  -g dataset/genes.gtf -f dataset/genes.fa
```





