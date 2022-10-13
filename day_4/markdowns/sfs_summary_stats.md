Site Frequency Spectrum & Summary Statistics
============================================

<br>

The following exercises will teach you how to estimate the site frequency spectrum (SFS) and population 
genetic summary statistics that are based on the SFS from low coverage whole-genome sequencing data using ANGSD.

### Initial Preparation

For today's exercises please generate the following directories

```bash
mkdir ~/day4
mkdir ~/day4/Results
mkdir ~/day4/Data
```

and set the following environment variables


```bash
DIR=/home/ubuntu/Share/data
DATA=$DIR/BAMS
REF=$DIR/Ref.fa
ANC=$DIR/outgrp_ref.fa
RESDIR=~/day4/Results
DATDIR=~/day4/Data
ANGSD=/home/ubuntu/Share/angsd/angsd
REALSFS=/home/ubuntu/Share/angsd/misc/realSFS
```

### Data

We will continue to work with the low-coverage NGS data for 60 Atlantic silversides from the following four populations:

<img src="https://github.com/nt246/physalia-lcwgs/blob/main/day_3/img/Silverside_Sample_Map.png?raw=true" height="400">

These populations have been previously studied in [Therkildsen et al. 2019](https://science.sciencemag.org/content/365/6452/487) and [Wilder et al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.189), and cover the entire range of the species.

Our NGS data are in BAM format and span a 2Mb region on chromosome 24. The interesting aspect about chromosome 24 is that it harbours a large polymorphic inversion that differs in its frequency across populations. 
The test dataset spans one breakpoint of this inversion (1Mb up and downstream).



### The Site Frequency Spectrum (SFS)

The SFS characterizes the distribution of allele frequencies in a population. That is, it records the proportion of sites in different allele 
frequency categories. This summary of allele frequencies is useful for assessing data quality, inferring demography, and detecting selection.
The "unfolded" SFS characterizes the frequency of derived alleles and requires some way to "polarize" alleles, i.e. decide which allelic state is 
ancestral and which is derived. It's common to use an outgroup sequence for this. The SFS can also be folded, in which case it characterizes the 
distribution of minor allele frequencies in a population. The folded SFS considers allele frequency classes of 1/2N to 0.5 (where N is the diploid sample size). 
To fold the spectrum sites with derived allele frequencies of (2N-1)/2N are in the same class as 1/2N sites, (2N-2)/2N sites are in the same class as a 2/2N 
site, (2N-3)/2N sites are in the same class as 3/2N sites, up to a class of 0.5 allele frequency (the highest frequency a minor allele can 
take by definition).

We will use ANGSD to estimate the SFS using the methods described [here](http://www.ncbi.nlm.nih.gov/pubmed/22911679).
Details on the implementation can be found [here](http://popgen.dk/angsd/index.php/RealSFSmethod). The main Wiki page describing the SFS and multidimensional 
SFS is [here](http://popgen.dk/angsd/index.php/SFS_Estimation).

The general workflow is pictured as

<img src="https://github.com/nt246/physalia-lcwgs/blob/main/day_3/img/Silverside_Sample_Map.png?raw=true" height="400">

We will estimate the unfolded SFS for the PANY population. To do this first we need to estimate the likelihood of sampling **k** derived alleles for k=0, k=1, k=2, ...,k=2N at 
every site. This is accomplished using `-doSaf`. We will start with BAMs as input and as usual we will use `-GL` to calculate genotype likelihoods from which we can start estimating the allele frequency likelihoods.

`-doSaf`

```bash
...
	-doSaf		0
	   1: SAF calculation integrating over possible minor alleles
	   2: SAF calculation incorporating inbreeding
	   3: Calculate genotype probabilities using SAF (DEPRECATED; use -doPost 3)
	   4: SAF calculation from genotype posteriors (input is beagle text format)
	   5: SAF calculation conditioning on minor allele from -doMajorMinor
	 -underFlowProtect	0
	 -anc		(null)	(ancestral fasta)
	 -noTrans	0	(remove transitions)
	 -pest		(null)	(prior SFS)
	 -isHap		0	(samples are haploid; works with -doSaf 1 or 5)
	 -scoreTol	1.0e-09	(tolerance for score-limited algorithm)
	 -doPost	0	(doPost 3, used for accessing SAF based variables)

NB: If -pest is supplied in addition to -doSaf then the output will be posterior probabilities of the sample allele frequency for each site
NB: Increasing -scoreTol will trade accuracy for reduced computation time and storage
```

Let's calculate the allele frequency likelihoods at all sites. You will typically use `-doSaf 1` unless your population is highly inbred 
in which case `-doSaf 2` is more accurate.

```bash
$ANGSD -b $DIR/PANY_bams.txt -ref $REF -anc $ANC -out $RESDIR/PANY \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 60 -doCounts 1 \
   -GL 1 -doSaf 1
```

You can have a look at the results using `realSFS` (note that realSFS is querying the index to make a human-readable 
output that combines the .saf.gz and .saf.pos.gz files):

```bash
$REALSFS print $RESDIR/PANY.saf.idx | less -S
```

The columns are (1) chromosome, (2) position, (3) likelihood of 0 derived alleles, (4) likelihood of 1 derived allele, 
(5), likelihood of 2 derived alleles, ..., (2N+2) likelihood of having 2N derived alleles. Note that the likelihoods 
are scaled to the highest likelihood and log transformed (such that the most likely allele frequency will have a 
value of 0).

*QUESTION*
Can you tell which are the first two sites that look to be variable?

<details>

<summary> Click for answer </summary>

Mme_chr24:2558528-4558528	48
Mme_chr24:2558528-4558528	61

</details>

*QUESTION*

What is the most likely number of derived alleles at position Mme_chr24:2558528-4558528 61?

<details>

<summary> Click for answer </details>

3 derived alleles.

</details>

Help information for the `realSFS` program can be seen by calling `realSFS` without arguments.

```bash
Usage:
	realSFS [options] deme1.saf.idx [deme2.saf.idx ...]	(calculates [multi-]SFS)
	realSFS subcommand					(displays usage for subcommands)
Options:
	-r chrom[:start-end]	(use only sites in specified region)
	-sites ?		(use only these sites, see 'Notes')
	-anc ?			(??, see 'Notes')
	-ref ?			(??, see 'Notes')
	-cores 1		(number of threads)
	-tole 1e-10		(convergence tolerance for EM)
	-maxiter 100		(maximum number of EM iterations)
	-bootstrap 0		(number of bootstrap replicates)
	-resample_chr 0		(0 = bootstrap sites; 1 = bootstrap chromosomes/contigs)
	-emaccl 1		(0 = regular EM; 1 = accelerated EM)
	-fold 0			(0 = unfolded SFS; 1 = folded SFS)
	-nsites 0		(number of sites to use in calculation; 0 = all sites)
	-seed -1		(random seed for start values; -1 = machine noise)
Subcommands:
	bins			(print SFS bins corresponding to flattened output)
	cat			(concatenate two SAF files)
	check			(checks that positions are ordered correctly)
	dadi			(call SFS bin for each site via empirical Bayes)
	fst			(index and calculate per-site Fst)
	print			(print SAF in various formats)
	print_header		(print SAF index information)
	saf2theta		(create inputs for theta calculation in bands)
	text2saf		(convert text saffile to binary safv3 file)
	winsfs			(Window SFS might be more stabile)
Examples:
	#one-dimensional SFS
	./realSFS deme.saf.idx

	#one-dimensional SFS for all of chromosome 22
	./realSFS deme.saf.idx -r chr22

	#two-dimensional SFS for all of chromosome 22
	./realSFS deme1.saf.idx deme2.saf.idx -r chr22

	#estimate the SFS for the first 500Mb (including multiple chromosomes)
	./realSFS deme.saf.idx -nsites 500000000

	#estimate the SFS for a specified region (e.g. around a gene)
	./realSFS deme.saf.idx -r chr2:135000000-140000000

	#generate 100 bootstrap replicates of SFS by resampling contigs
	./realSFS deme.saf.idx -bootstrap 100 -resample_chr 1
Notes:
	Output is the maximum likelihood estimate of sites for each SFS bin, as a flattened array. If multiple index files are supplied, the joint SFS is calculated. To see the SFS bins associated with each value in the array, use the 'bins' subcommand (especially useful for folded multidimensional SFS). Bootstrap replicates are output one/line after the MLE.
```

Now we will use the allele frequency likelihoods calculated for every site contained in the .saf file to obtain 
a maximum likelihood estimate of the SFS with the `realSFS` program.

```bash
$REALSFS $RESDIR/PANY.saf.idx > $RESDIR/PANY.sfs
```

Take a look at the output:

```bash
cat $RESDIR/PANY.sfs
```
You should observe the following vector of 2N+1 values:

```bash
1287585.159170 5441.014450 1893.638927 1024.182258 1477.927761 1192.652825 353.801352 166.179066 199.498879 736.762268 775.654171 132.182867 36.552615 55.659970 162.170177 327.757461 287.965161 106.914226 46.408186 75.780528 238.957573 419.965002 305.941344 92.168706 38.235800 89.032890 167.489915 426.555830 265.498645 236.967396 47411.324580
```


These are the *expected* number of sites in the PANY sample with 0 (value 1), 1 (value 2), 3 (value 3), ..., 2N (value 2N+1) derived alleles.

Plotting the SFS (TO DO)
```
/usr/bin/R

sfs<-scan("smallFolded.sfs")

pdf("~/day4/Results/PANY_1D_SFS.pdf")
plot = barplot(sfs[-1])
dev.off() 
```




What if we did not have an outgroup to polarize alleles with? In this case we could calculate the folded SFS. Let's try it.

First we'll calculate the allele frequency likelihoods with `-doSaf 1`, which requires us to supply `-anc` to which we will just pass 
the reference FASTA in the case of folding (remember, we are pretending that we don't have a reliable ancestral fasta to polarize with):


```bash
$ANGSD -b $DIR/PANY_bams.txt -ref $REF -anc $REF -out $RESDIR/PANY_fold \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 60 -doCounts 1 \
   -GL 1 -doSaf 1
```

Now calculate the SFS using this new .saf file as before, except this time tell `realSFS` to fold the spectrum with `-fold 1`

```bash
$REALSFS -fold 1 $RESDIR/PANY_fold.saf.idx > $RESDIR/PANY_fold.sfs
```

Let's look at the ML estimate of the folded SFS:

```bash
cat $RESDIR/PANY_fold.sfs
```

Now your vector of expected counts should look like this:
1403715.239733 6987.037568 2279.638485 2180.218338 1813.594151 1117.329380 697.428766 652.617721 797.370326 889.151374 744.767130 541.466942 478.770881 458.365587 405.231595 187.772023 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000

The folded SFS is represented by the first 1 through N values. These are the expected number of sites with minor allele counts of 0, 1, 2, ..., N. `realSFS` still prints 2N+1 values but 
the N+1 to 2N+1 values should all be `0.000000` and can be ignored when using `-fold 1` because these categories don't really exist in the folded SFS.

### Allele frequency posterior probabilities

Since you were able to calculate allele frequency likelihoods with `-doSaf` you may be wondering if you can use these to calculate the posterior probabilities of the possible 
allele frequencies at each site. You can! In order to do this you can use the site frequency spectrum calculated with `realSFS` as a prior on the allele frequencies since it tells 
you the probability of randomly drawing a site with a given number of derived (or minor if folded) alleles in the population.

Let's calculate derived allele frequency posterior probabilities for the PANY population. This is achieved by running the per site allele frequency likelihood calculation (`-doSaf`) while supplying 
the SFS are a prior with '-pest':

```bash
$ANGSD -b $DIR/PANY_bams.txt -ref $REF -anc $ANC -out $RESDIR/PANY_post \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 60 -doCounts 1 \
   -GL 1 -doSaf 1 -pest $RESDIR/PANY.sfs
```bash

Have a look at the output:

```bash
$REALSFS print $RESDIR/PANY_post.saf.idx | less -S
```
The columns are (1) chromosome, (2) position, and then the posterior probability of having 0, 1, 2, 3, ..., 2N derived alleles in each of the subsequent columns. The posterior 
probablities are log transformed.

*QUESTION*

What is the most probable number of derived alleles at Mme_chr24:2558528-4558528 61?

<details>

<summary> Click for answer </details>

0 derived alleles.

</details>

*QUESTION*

What is the posterior probability that position Mme_chr24:2558528-4558528 61 has 1 derived allele (i.e. is a singleton)?

<details>

<summary> Click for answer </summary>

The second value in the probabillity vector for Mme_chr24:2558528-4558528 61 is the probability of 1 derived allele in log space. 
This value is -1.844222. So the probability that this site is a singleton is exp(-1.844222) = *0.1581483*.

</details>

*QUESTION*

What is the probability that the following sites are variable?

Mme_chr24:2558528-4558528	39
Mme_chr24:2558528-4558528	48
Mme_chr24:2558528-4558528	61

<details>

<summary> Click for answer </summary>

The probabiliyt that a site is variable is given by, P(variable) = 1 - P(0 derived alleles) + P(2N derived alleles). This is because sites with either 0 or 2N derived 
alleles are fixed.

P(Mme_chr24:2558528-4558528 39 is variable) = 1 - exp(-0.005389) + exp(-Inf) = 0.005374505
P(Mme_chr24:2558528-4558528 48 is variable) = 1 - exp(-Inf) + exp(-Inf) = 1
P(Mme_chr24:2558528-4558528 61 is variable) = 1 - exp(-0.638232) + exp(-Inf) = 0.4717745

Alternatively, you could take the sum over P(**x** derived alleles) for **x**=1 to **x**=2N-1. This would give the same answers.

</details>

*BONUS QUESTION*

What would the SFS calculated for a single diploid individual tell you?

<details>

<summary> Click for answer </summary>

The possible allele frequency categories for a single diploid individual are 0, 1, and 2 derived alleles for the unfolded SFS  
or just 0 and 1 minor alleles for the folded SFS. A site with one derived (or minor) allele in a single individual represents a heterozygous site, 
so the SFS for a single individual provides the heterozygosity of that individual.

</details>

### Diversity and Neutrality Statistics

Many population genetic summary statistics can be calculated directly from the SFS.

./misc/realSFS saf2theta out.saf.idx -outname out -sfs out.sfs (TODO)

### Multidimensional SFS and Fst

(TODO)
