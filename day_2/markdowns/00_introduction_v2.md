
In this session you will learn how to:
* Calculate genotype likelihoods and call genotypes in a probabilistic framework
* Estimate allele frequencies
* Call SNPs in the face of high genotyping uncertainty

We are using the program [ANGSD](http://popgen.dk/angsd/index.php/ANGSD) (Analysis of Next Generation Sequencing Data).
More information about the methods implemented in ANGSD can be found in its associated [paper](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

*ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input 
types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. 
This is especially useful for low and medium depth data.* (quoted from the ANGSD website).

Please make sure to follow the preparatory instructions below before continuing with the exercises. 
You need to set the path to the software and various data that will be used.
You also need to create two folders in your home user directory; one for your results and one for your intermediate data.

Create a folder for this session in your home directory and enter it.
```
mkdir ~/day2
cd ~/day2

```
You should now be in `~/day2`. You can confirm this with `pwd`, which will print your location in the filesystem.

Now create two folders in your home user directory, one for your results and one for your intermediate data.
```
mkdir ~/day2/Results
mkdir ~/day2/Data

```
Set environment variables to point to these directories, files that we'll use, and software.
```
RESDIR=~/day2/Results
DATDIR=~/day2/Data
DIR=/home/ubuntu/Share/physalia-lcwgs/data
DATA=$DIR/BAMS
REF=$DIR/Ref.fa
ANC=$DIR/outgrp_ref.fa
angsd=/home/ubuntu/angsd/angsd

```

The **workflow** is roughly divided into four steps:

[1. Data filtering options and I/O](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/01_filtering.md)

[2. Genotype likelihoods](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/02_likelihoods.md)

[3. Allele frequency estimation & SNP calling](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/03_allele_frequencies.md)

[4. Genotype calling](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/04_genotype.md)

You are now going to learn how to build your first ANGSD command for data processing and filtering.

[click here](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/01_filtering.md) to move to the next session.

-----------------------------------------------


