# Populaiton structure analysis

In this session you will learn how to do:

* Principal Components Analysis
* Admixture analysis

using low-coverage whole genome data.

We have low-coverage NGS data for 60 Atlantic silversides from four populations, spanning a 2Mb region on chromosome 24. These populations have been previously studied in [Therkildsen et al. 2019](https://science.sciencemag.org/content/365/6452/487) and [Wilder et al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.189) using a trancsriptome as reference, and cover the entire distribution range of Atlantic silverside. 
The interesting aspect about chromosome 24 is that it harbours a large polymorphic inversion, and our test dataset spans one breakpoint of this inversion. 

First, can we determine the neutral population structure of Atlantic silverside along its distribution range. 
Second, can use principal components analysis and ancestry analyses to determine which populations contain the alternative inversion karyotype? And are there potentially any heterozygous individuals for the inversion?

For this practical, we will be using [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data), [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd).


Please make sure to follow these preparatory instructions below before running these examples. 
Briefly, you need to set the path to the software and various data that will be used.
Also, you will have to create two folders on your working directory, one for your results and one for your intermediate data.

```
NGS=/programs/angsd0.930/
DIR=/workdir/arne/physalia_lcwgs_data/data_practicals/
DATA=$DIR/BAMS/
REF=$DIR/Ref.fa

mkdir Results
mkdir Data
```

--------------------------------------------------


## Principal components analysis (PCA)

To perform PCA with low-coverage NGS data, we have to infer the genetic covariance matrix, which can be estimated in different ways.
Here we will estimate the covariance matrix using single-read sampling in [ANGSD](http://www.popgen.dk/angsd/index.php/PCA_MDS) but will also discuss the algorithm implemented in [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsd)

In general it is good practice to perform LD pruning or at least thinning to reduce the impact of large LD clusters (e.g. inversions) on the inferred population structure. In this practical, we will perform the population structure analysis using an LD-pruned SNP dataset that you prepared earlier today and also the full SNP dataset. The full SNP dataset will help us to understand how the inversion karyotype is distributed in our dataset and will highlight the impact of large LD clusters on the inferred structure. 

For the PCA method you should use only inferred variant sites (SNPs). SNPs can be inferred based on genotype likelihoods (see SNP_calling day 2) at the same time as inferring the covariance matrix (`-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1`) by providing the `-SNP_pval` option. SNPs should be restricted to more common variants with minor allele frequencies of at least 5% using the `-minMAF 0.05` option. However, this way, you will focus on all SNPs and not only an LD-pruned subset.

```
# $NGS/angsd/angsd -b ALL_bams.txt -anc $REF -out Results/MME_ANGSD_PCA \
#	-minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
#	-GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
#	-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -P 4
```

Alternatively, you can provide a list of variant sites using the `-sites` options. Here we will use a list of LD-pruned variant sites created earlier today. 



### PCA with LD-pruned SNP dataset

Using the list of LD-pruned variant sites and the code shown above, we can estimate a covariance matrix using single-read sampling for LD-pruned SNPs using ANGSD. 

```
$NGS/angsd/angsd -b ALL_bams.txt -anc $REF -out Results/MME_ANGSD_PCA_LDpruned \
	-minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
	-GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
	-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -sites LDpruned_snps.list
```

At the same time, we will also output the genotype likelihoods for these variant sites in beagle likelihood file format (beagle.gz), which will be used as input for estimating the covariance matrix using PCAngsd (below) 

```
library(tidyverse) #load the tidyverse package for formatting and plotting

#Load the covariance matrix
cov <- as.matrix(read.table("/workdir/arne/physalia_lcwgs_data/data_practicals/Results/MME_ANGSD_PCA_LDpruned.covMat", header = F))

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors <- (mme.pca$vectors) #extract eigenvectors 
pca.vectors <- as_tibble(cbind(pop, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors) #check all columns and convert if necessary (automatic)

#plot PC1 vs PC2 using ggplot
pca = ggplot(data = df, aes(x=V2, y=V3, fill = pop, colour = pop)) +
  geom_point(size = 5, shape = 21) +
  xlab("Principal component 1") +
  ylab("Principal component 2")

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_plot.pdf", plot = pca)

```
You can run the Rscript containing this code to plot your results:
```
Rscript ./plotPCA.R  
```


### Optional: PCA with all SNPs

If you are interested, you can compare the PCA with LD-pruned SNPs to one that was performed for all SNPs without LD-pruning.

What is the difference in the PCA pattern? 

1. Covariance matrix and PCA with ANGSD 

```
$NGS/angsd/angsd -b ALL_bams.txt -anc $REF -out Results/MME_ANGSD_PCA \
	-minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
	-GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
	-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -P 4
```

Now we can use the covariance matrix to perform a principal components analysis using the `eigen` function in R:

```
library(tidyverse) #load the tidyverse package for formatting and plotting

#Load the covariance matrix
cov <- as.matrix(read.table("/workdir/arne/physalia_lcwgs_data/data_practicals/Results/MME_ANGSD_PCA.covMat", header = F))

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors <- (mme.pca$vectors) #extract eigenvectors 
pca.vectors <- as_tibble(cbind(pop, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors) #check all columns and convert if necessary (automatic)

#plot PC1 vs PC2 using ggplot
pca = ggplot(data = df, aes(x=V2, y=V3, fill = pop, colour = pop)) +
  geom_point(size = 5, shape = 21) +
  xlab("Principal component 1") +
  ylab("Principal component 2")

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_plot.pdf", plot = pca)

```
You can run the Rscript containing this code to plot your results:
```
Rscript ./plotPCA.R  
```

Comparing the PCA to the expected pattern based on the whole genome (see lecture) we see that these patterns differ drastically.
The difference in population structure can be explained by the low number of SNPs we have in our test dataset. 


2. Optional: Covariance matrix estimation with PCAngsd (based on all SNPs)

A covariance matrix can also be inferred from genotype likelihoods using PCAangsd. 
PCAngsd takes as input genotype likelihoods in beagle format. 
We have already created the beagle input file in the previous step using the `-doGLF 2` option while we created the covariance matrix using the read re-sampling.  


GO THROUGH POTENTIAL SETTINGS:

```
python3 /programs/pcangsd-0.98/pcangsd.py -beagle Results/MME_ANGSD_PCA.beagle.gz -o Results/covmatrix
```

From PCAngsd 0.98 the output is saved in numpy format but this can easily loaded into R (see script below).
We can perform the principal components analysis and plot PC1 vs PC2 the same way we did before.

```
library(tidyverse) #load the tidyverse package for formatting and plotting

#Load the covariance matrix
cov = as.matrix(RcppCNPy::npyLoad("/Users/arnejacobs/Dropbox/covmatrix.cov.npy"))

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

mme.pca <- eigen(cov) #perform the pca using the eigen function. 

eigenvectors <- (mme.pca$vectors) #extract eigenvectors 
pca.vectors <- as_tibble(cbind(pop, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors) #check all columns and convert if necessary (automatic)

#plot PC1 vs PC2 using ggplot
pca = ggplot(data = df, aes(x=V2, y=V3, fill = pop, colour = pop)) +
  geom_point(size = 5, shape = 21) +
  xlab("Principal component 1") +
  ylab("Principal component 2")

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_pcangsd_plot.pdf", plot = pca)
```

You can run the Rscript containing this code to plot your results:
```
Rscript ./plot_pca_pcangsd.R  

```


--------------------------------------------------

## Admixture analysis

In some cases we also want to infer admixture proportions. Similar to PCA, there are different ways of inferring admixture proportions from genotype likelihoods.
We will go in more detail through [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and will present a way of inferring admixture proportions with PCAngsd.  

1. ngsAdmix

A genotype likelihood file in beagle format (same as for PCAngsd) is specified as input using the `-likes` option. 
In addition, there are a range of parameters that can be adjusted. Here we only set the number of ancestry clusters using the `-K` option.
In reality, it is advisable to compare different numbers of ancestry clusters by changing the `-K` option. 

Also, if the analysis is not converging, one might want to increase the maximum number of EM iterations using the `-maxiter` option. 

Here we will run ngsAdmix with one value of K. 
```
/programs/NGSadmix/NGSadmix -likes Results/MME_ANGSD_PCA.beagle.gz -K 2 -o Results/MME_ngsAdmix_K2_out
```

Plot admixture proportions in R e.g. following this Rscript:

```
Rscript plot_ngsAdmix.R
```


If you want to, you can try changing the value for `-K` and compare the infferred likelihoods and results. 
```
/programs/NGSadmix/NGSadmix -likes Results/MME_ANGSD_PCA.beagle.gz -K 3 -o Results/MME_ngsAdmix_K3_out
```


2. PCAngsd

In addition to ngsAdmix, PCAngsd can also infer admixture proportions. 
The command is similar to the one used for the PCA with the addition of to a `-admix` option. 
The input file for ngsAdmix and PCAngsd is the same, making comparisons relatively easy. 

Other than ngsAdmix, PCAngsd can automatically infer the most likely number of ancestry cluster.
However, one can also set the number of clusters using the `-admix_K` option. 


```
python3 /programs/pcangsd-0.98/pcangsd.py -beagle Results/MME_ANGSD_PCA.beagle.gz -admix -admix_K 2 -o Results/MME_PCAngsd_K2_out

```




















