# Populaiton structure analysis

In this session you will learn how to do:

* Principal Components Analysis
* Admixture analysis

using low-coverage whole genome data.

We have low-coverage NGS data for 60 Atlantic silversides from four populations, spanning a 2Mb region on chromosome 24. These populations have been previously studied in [Therkildsen et al. 2019](https://science.sciencemag.org/content/365/6452/487) and [Wilder et al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.189), and cover the entire distribution range of Atlantic silverside. 
The interesting aspect about chromosome 24 is that it harbours a large polymorphic inversion that differs in its frequency across populations, and our test dataset spans one breakpoint of this inversion (1Mb up and downstream). Therefore, we might expect that differ parts of our dataset contain different evolutionary signals.

In this tutorial, can want to use these data to determine the neutral population structure of Atlantic silverside along its distribution range. 
Second, we can use principal components analysis and ancestry analyses to determine which populations contain the alternative inversion karyotype? And are there potentially any heterozygous individuals for the inversion?

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

In general it is good practice to perform LD pruning or at least thinning to reduce the impact of non-independent SNP clusters, e.g. large LD clusters in inversions, on the inferred population structure. In this practical, we will perform the population structure analysis using an LD-pruned SNP dataset that you prepared earlier today and also the full SNP dataset. The full SNP dataset will help us to understand how the inversion karyotype is distributed in our dataset and will highlight the impact of large LD clusters on the inferred structure. 

For the PCA method you should use only inferred variant sites (SNPs). SNPs can be inferred based on genotype likelihoods (see SNP_calling day 2) at the same time as inferring the covariance matrix (`-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1`) by providing the `-SNP_pval` option. SNPs should be restricted to more common variants with minor allele frequencies of at least 5% using the `-minMAF 0.05` option. However, this way, you will focus on all SNPs and not only an LD-pruned subset.

```
# $NGS/angsd/angsd -b ALL_bams.txt -anc $REF -out Results/MME_ANGSD_PCA \
#	-minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
#	-GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
#	-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -P 4
```

Alternatively, you can provide a list of variant sites using the `-sites` options. Here we will use a list of LD-pruned variant sites created earlier today. 



### PCA with LD-pruned SNP dataset


**1. Estimating the covariance matrix using a single read sampling approach:**

Using the list of LD-pruned variant sites and the code shown above, we can estimate a covariance matrix using single-read sampling for LD-pruned SNPs using ANGSD. 

```
$NGS/angsd/angsd -b ALL_bams.txt -anc $REF -out Results/MME_ANGSD_PCA_LDpruned \
	-minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
	-GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
	-doCounts 1 -doCov 1 -sites LDpruned_snps.list
```

At the same time, we will also output the genotype likelihoods for these variant sites in beagle likelihood file format (beagle.gz), which will be used as input for estimating the covariance matrix using PCAngsd (below) 


The actual principal components analysis is then performed in R using the `eigen` function. For that, we have to load the covariance matrix into R, then we can optionally provide population assignments for each individual (rows in same order as input bam file list `-b ALL_bams.txt`).

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
```

We can then extract the eigenvectors from the pca object and format them into a dataframe for plotting, e.g. using `ggplot()`.
```
eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, eigenvectors)) #combine with our population assignments
df = type_convert(pca.vectors) #check all columns and convert if necessary (automatic)

#plot PC1 vs PC2 using ggplot
pca = ggplot(data = df, aes(x=V2, y=V3, fill = pop, colour = pop)) +
  geom_point(size = 5, shape = 21) +
  xlab("Principal component 1") +
  ylab("Principal component 2")

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_plot.pdf", plot = pca)
```

Additionally, we can extract the eigenvalues for each eigenvector, and can then estimate the variance explained for each eigenvector (e.g. here for PC1 to PC4): 
```
pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
```

You can run this Rscript to estimate and plot the principal components:
```
Rscript ./plotPCA.R  
```


**2. Alternative approach: Covariance matrix estimation with PCAngsd**

The covariance matrix can also be inferred from genotype likelihoods using PCAangsd. PCAngsd takes as input genotype likelihoods in beagle format, which we generated in the step before using the `-doGLF 2` option.

PCAngsd provides a multitude of different settings, described [here](http://www.popgen.dk/software/index.php/PCAngsd). We won't change any of the settings here and only use the default settings, which are sufficient in most cases. 

We provide the path to the input file using the `-beagle` option, which also tells PCAngsd that we are working with a beagle file. The output path and output name is provided using the `-o` option. 

```
python3 /programs/pcangsd-0.98/pcangsd.py -beagle Results/MME_ANGSD_PCA.beagle.gz -o Results/covmatrix
```

From PCAngsd 0.98, the output is saved in numpy format but this can easily loaded into R using the `RcppCNPy::npyLoad` R function (see script below).

We can perform the principal components analysis and plot PC1 vs PC2 the same way we did before, except that the function for importing the covariance matrix has to be adjusted.

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

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_pcangsd_plot.pdf", plot = pca)
```

You can run the Rscript containing this code to plot your results:
```
Rscript ./plot_pca_pcangsd.R  

```



### Optional: Compare the results to a PCA based on all SNPs

While it is best practice to perform a PCA based on an LD-pruned SNP dataset, PCAs are often performed based on full SNP datasets. 
If you are interested and have time, you could compare the PCA for LD-pruned SNPs to one that was performed for all SNPs in the dataset without LD-pruning.

Q: What is the difference in the PCA pattern? 

For this, we will only estimate the covariance matrix using single read sampling in ANGSD:

```
$NGS/angsd/angsd -b ALL_bams.txt -anc $REF -out Results/MME_ANGSD_PCA \
	-minMapQ 20 -minQ 20 -doMaf 1 -minMaf 0.05 -SNP_pval 2e-6 \
	-GL 1 -doGlf 2 -doMajorMinor 1 -doPost 1 \
	-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -P 4
```

Again, we perform the principal components analysis using the `eigen` function in R:
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
You can run the plot your results using this script:
```
Rscript ./plotPCA.R  
```

Q: How do population relationships differ between the on LD-pruned PCA and the full-dataset PCA? And can you guess which individuals carry which inversion karyotype? 





--------------------------------------------------

## Admixture analysis

In some cases we also want to infer genome-wide admixture proportions for each individuals. Similar to PCA, there are different ways of inferring admixture proportions from genotype likelihoods. Here, we will use [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix) and will also present a way of inferring admixture proportions with PCAngsd.
Similar to the PCA, we want to use an LD-pruned SNP dataset for this analysis to reduce the impact of non-independent SNPs on the ancestry inference. 

### ngsAdmix

ngsAdmix uses a genotype likelihood file in beagle format (same as for PCAngsd) as input, which is specified using the `-likes` option. 
In addition, there are a range of parameters that can be adjusted. Here we only set the number of ancestry clusters using the `-K` option to K=2.
In reality, it is advisable to compare different numbers of ancestry clusters by iterating over different values of K. 

```
/programs/NGSadmix/NGSadmix -likes Results/MME_ANGSD_PCA.beagle.gz -K 2 -o Results/MME_ngsAdmix_K2_out
```

In case the analysis is not converging, one can also increase the maximum number of EM iterations using the `-maxiter` option. 

ngsAdmix produces three different outputs files:

* A run log containing the log likelihood of the admixture estimates: .log file
* A file containing the allele frequencies of all ancestral populations (in this case two ancestral clusters): .fopt file
* A file containing the admixture proportions for each individual: .qopt file

We are mostly interested in the admixture proportions and the log likelihood of the estimates. The log likelihoods can be compared between runs with different values of K to select the most likely number of ancestral clusters. 

You can plot admixture proportions in R e.g. following this Rscript:
```
Rscript plot_ngsAdmix.R
```
Description of the script:
```
library(tidyverse) #load the tidyverse package for formatting and plotting

#Load the covariance matrix
cov = read_table("/workdir/arne/physalia_lcwgs_data/data_practicals/Results/MME_ngsAdmix_K2_out.qopt", col_names = F)

#We will also add a column with population assingments
pop <- c("JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA","JIGA"
         ,"PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY","PANY"
         ,"MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS","MBNS"
         ,"MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU","MAQU")

cov.id = as.data.frame(cbind(pop, cov))
names(cov.id) = c("pop","q1","q2")
barplot(t(as.matrix(subset(cov.id, select=q1:q2))), col=c("firebrick","royalblue"), border=NA)

#Save plot as pdf
ggsave(filename = "/workdir/arne/physalia_lcwgs_data/data_practicals/Results/pca_pcangsd_plot.pdf", plot = pca)
```

If you want, you can try changing the value for `-K` and compare the infferred log ikelihoods and admixture proportions. Which value of K has the stronger support?
```
/programs/NGSadmix/NGSadmix -likes Results/MME_ANGSD_PCA.beagle.gz -K 3 -o Results/MME_ngsAdmix_K3_out
```


### Alternative approach: PCAngsd

In addition to ngsAdmix, PCAngsd can also be used to infer admixture proportions. 
The command is similar to the one used for the PCA with the addition of to a `-admix` option. 
The input file for ngsAdmix and PCAngsd is the same, making comparisons relatively easy. 

Other than ngsAdmix, PCAngsd can automatically infer the most likely number of ancestry cluster.
However, one can also set the number of clusters using the `-admix_K` option. 

```
python3 /programs/pcangsd-0.98/pcangsd.py -beagle Results/MME_ANGSD_PCA.beagle.gz -admix -admix_K 2 -o Results/MME_PCAngsd_K2_out

```





















