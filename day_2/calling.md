In this session you will learn how to do:
* genotype calling
* allele frequency estimation
* variant (or SNP) calling

We are using the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data).
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).

According to its website *ANGSD is a software for analyzing next generation sequencing data. The software can handle a number of different input types from mapped reads to imputed genotype probabilities. Most methods take genotype uncertainty into account instead of basing the analysis on called genotypes. This is especially useful for low and medium depth data.*

Please make sure to follow these preparatory instructions below before running these examples.
Briefly, you need to set the path to the software and various data that will be used.
Also, you will have to create two folders on your working directory, one for your results and one for your intermediate data.
```

NGS=/ricco/data/matteo/Software/ngsTools

DIR=/home/matteo/Copenhagen
DATA=/ricco/data/matteo/Data
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz

mkdir Results
mkdir Data
```

--------------------------------------------------

### Workflow

The workflow for this practical looks like this

![stages](calling_files/stages.png)

which seems daunting!
However, that's not the case and we will go through each step to understand each one of them.

The workflow is roughty divided into four steps:

0. Data filtering and I/O
1. Genotype likelihoods
2. Genotype calling
3. SNP calling

---------------------------------------------






