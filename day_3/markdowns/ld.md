Linkage disequilibrium estimation
================

- <a href="#define-paths-to-the-project-directory-and-programs"
  id="toc-define-paths-to-the-project-directory-and-programs">Define paths
  to the project directory and programs</a>
  - <a href="#set-the-project-directory-as-a-variable-named-basedir"
    id="toc-set-the-project-directory-as-a-variable-named-basedir">Set the
    project directory as a variable named <code>BASEDIR</code></a>
  - <a href="#specify-the-path-to-required-programs-as-variables"
    id="toc-specify-the-path-to-required-programs-as-variables">Specify the
    path to required programs as variables</a>
- <a href="#estimate-ld" id="toc-estimate-ld">Estimate LD</a>
  - <a href="#prepare-the-input-files"
    id="toc-prepare-the-input-files">Prepare the input files</a>
  - <a href="#run-ngsld" id="toc-run-ngsld">Run ngsLD</a>
- <a href="#visualize-ld-blocks" id="toc-visualize-ld-blocks">Visualize LD
  blocks</a>
- <a href="#ld-pruning" id="toc-ld-pruning">LD pruning</a>
  - <a href="#run-ld-pruning" id="toc-run-ld-pruning">Run LD pruning</a>
  - <a href="#generate-an-ld-pruned-snp-list"
    id="toc-generate-an-ld-pruned-snp-list">Generate an LD-pruned SNP
    list</a>
- <a href="#practical-considerations"
  id="toc-practical-considerations">Practical considerations</a>

<br> <br>

The estimation of linkage disequilibrium (LD) has important
applications, e.g. for inference of population size, demographic
history, selection, and for the discovery of structural variants. In
addition, since many downstream analyses make assumptions about the
independence of genomic loci, LD estimates are essential for trimming
the list of loci to be included in these analyses (LD pruning). In this
session, you will learn to estimate linkage disequilibrium using the
program [ngsLD](https://github.com/fgvieira/ngsLD), which employs a
maximum likelihood approach to account for genotype uncertainty in
low-coverage whole genome sequencing data. You will then visualize the
LD pattern and generate a list of LD-pruned SNPs for downstream
analyses.

We will continue working with low-coverage NGS data for 60 Atlantic
silversides from different four populations studied in [Therkildsen et
al. 2019](https://science.sciencemag.org/content/365/6452/487) and
[Wilder et
al. 2020](https://onlinelibrary.wiley.com/doi/10.1002/evl3.189). The
sequencing data from these individuals have been mapped to a 2Mb section
of chromosome 24. The entire genome is \~650 Mb, so this is just a small
snippet, but it’s an interesting one because it harbours a large
polymorphic inversion that varies substantially in frequency across the
species distribution range. Our test dataset spans one breakpoint of
this inversion (1Mb up and downstream). Therefore, we might expect that
LD patterns dramatically differ across the inversion breakpoint in our
dataset, as inversions can strongly suppress recombination.

<br>

## Define paths to the project directory and programs

We need to make sure the server knows where to find the programs we’ll
be running and our input and output directories. This will always need
to be specified every time we run our scripts in a new login session.

<br>

### Set the project directory as a variable named `BASEDIR`

We’ll copy today’s data files into your home directory with the
following commands

``` bash
## Copy the shared project directory to your home directory
cp -r ~/Share/physalia-lcwgs/day_3 ~/day3

## Define BASEDIR as your project directory
BASEDIR=~/day3/ # Note that no spaces are allowed! 

cd $BASEDIR
ls
```

<br>

### Specify the path to required programs as variables

When running these scripts on the AWS server, run the following:

``` bash
NGSLD=~/Share/ngsLD/
## launch the conda environment that contains the dependencies of ngsLD submodules
conda activate ngsLD
## add library paths
export LIBRARY_PATH=/home/ubuntu/src/conda/envs/ngsLD/lib/
export LD_LIBRARY_PATH=/home/ubuntu/src/conda/envs/ngsLD/lib/:$LD_LIBRARY_PATH
```

<br>

If you will be running these programs on a different system, you will
have to specify the paths to the different programs on that system (or
add them to your \$PATH).

<br>

## Estimate LD

### Prepare the input files

ngsLD requires two input files.

1.  `--geno FILE`: input file with genotypes, genotype likelihoods or
    genotype posterior probabilities. With low-coverage data, a genotype
    likelihood file is often used, in which each row is a SNP and each
    individual has three columns corresponding to the likelihood of the
    three genotypes (we only need to keep track of the likelihood of
    three genotypes (rather than all ten possible) when the major and
    minor allele has been inferred). Therefore, a `beagle` formatted
    genotype likelihood file generated from ANGSD (`-doGlf 2`) can be
    inputted into ngsLD after the header row and the first three columns
    (i.e. positions, major allele, minor allele) are removed. Here,
    because of time constraint, we will subsample a beagle file
    generated from a similar ANGSD run to what we explored yesterday
    (`MME_ANGSD_PCA.beagle.gz`) as the input to ngsLD by selecting one
    SNP in every 50 SNPs.

2.  `--pos FILE`: input file with site coordinates (one per line), where
    the 1st column stands for the chromosome/contig and the 2nd for the
    position (bp). One convenient way to generate this is by selecting
    the first two columns of the `mafs` file outputted by ANGSD, with
    the header removed. Again, for this exercise, we will downsample the
    `mafs` file that were generated in the same ANGSD run that generated
    our genotype likelihood beagle file (`MME_ANGSD_PCA.mafs.gz`).

``` bash
## Prepare a geno file by subsampling one SNP in every 50 SNPs in the beagle file
zcat $BASEDIR/angsd/MME_ANGSD_PCA.beagle.gz | awk 'NR % 50 == 0' | cut -f 4- | gzip  > $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.beagle.gz

## Prepare a pos file by subsampling one SNP in every 50 SNPs in the mafs filre
zcat $BASEDIR/angsd/MME_ANGSD_PCA.mafs.gz | cut -f 1,2 |  awk 'NR % 50 == 0' | sed 's/:/_/g'| gzip > $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.pos.gz
```

Don’t worry if you don’t understand the `sed` command. The `:` in the
chromosome name interferes with the code, so the `sed` command just
replaces the `:` with `_`.

You can use something like
`zcat $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.beagle.gz | cut -f 1,2,3 | head`
and `zcat $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.pos.gz | head` to
check the content of these input files.

<br>

### Run ngsLD

Important ngsLD parameters (default values shown in square brackets):

- `--probs`: specification of whether the input is genotype
  probabilities (likelihoods or posteriors).
- `--n_ind INT`: sample size (number of individuals).
- `--n_sites INT`: total number of sites.
- `--max_kb_dist DOUBLE`: maximum distance between SNPs (in Kb) to
  calculate LD. Set to 0 (zero) to disable filter. \[100\]
- `--max_snp_dist INT`: maximum distance between SNPs (in number of
  SNPs) to calculate LD. Set to 0 (zero) to disable filter. \[0\]
- `--n_threads INT`: number of threads to use. Using more threads will
  speed up the computation. \[1\]
- `--out FILE`: output file name. \[stdout\]

``` bash
$NGSLD/ngsLD \
--geno $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.beagle.gz \
--pos $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.pos.gz \
--probs \
--n_ind 60 \
--n_sites 1134 \
--max_kb_dist 0 \
--n_threads 1 \
--out $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld 
```

This step may take a few minutes to run. In the meantime, you can take a
look at the [ngsLD GitHub page](https://github.com/fgvieira/ngsLD).
Particularly, read [this
section](https://github.com/fgvieira/ngsLD#output) on output format, and
see if you can identify the different measures of LD in the output file
once it is generated. Note that the r<sup>2</sup> value from the EM
algorithm is the one that is going to be the most useful.

<br>

## Visualize LD blocks

We will use a script slightly modified from the original `LD_blocks.sh`
script provided by ngsLD to generate a plot of LD blocks in our data. It
takes three argument in the following order:

1.  chromosome / linkage group / scaffold name
2.  starting position
3.  ending position

``` bash
cd $BASEDIR/ngsld
cat $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld | bash $BASEDIR/scripts/LD_blocks.sh \
Mme_chr24_2558528-4558528 \
200000 \
1400000
```

To view the resulting plot, you can use `scp` to transfer the pdf file
to your local computer. Note that you will need to edit the .pem file
path and name, user name, IP address, and the destination path in the
following script, and run it on your local computer \[i.e. not in the
Terminal window you are logged into the server with\].

``` bash
## scp -i c1.pem user1@35.164.117.54:/home/user1/day3/ngsld/LD_blocks.r2.pdf ./
```

<br>

Examine the pdf file that the script outputs. If you can not get the
code to work, see the [plot that we have
generated](https://github.com/nt246/physalia-lcwgs/blob/main/day_3/ngsld/LD_blocks.r2.pdf).

Do you notice any interesting pattern in this plot of LD blocks? What do
you think is causing this pattern and why?

<br>

## LD pruning

### Run LD pruning

For many downstream analyses, independence among different SNPs is often
assumed, so it is important to generate a list of SNPs that are in low
LD with each other. To do this, we can use the `prune_ngsLD.py` script
provided in ngsLD. This script takes the LD estimation output (in our
case `MME_ANGSD_PCA.ld`) as its input. Some important parameters
include:

- `--input FILE`: Path to input file (i.e. `.ld` file generated from the
  LD estimation step)
- `--max_dist INT`: Maximum distance between two SNPs to assume they are
  potentially in strong LD
- `--min_weight FLOAT`: Minimum LD estimates between two SNPs to assume
  that they are in strong LD (and thus one of them will be pruned)
- `--output FILE`: Path to output file \[STDOUT\]

<br>

Run the following script **on the AWS server** to perform LD pruning
with our test dataset.

``` bash
## remove the header line from LD estimates (because the python script had issues with it)
tail -n +2 $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.ld > $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.no_header.ld
## run LD pruning script
python $NGSLD/scripts/prune_ngsLD.py \
--input $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled.no_header.ld \
--max_dist 2000000 \
--min_weight 0.5 \
--output $BASEDIR/ngsld/MME_ANGSD_PCA_subsampled_unlinked.id
```

Check the LD pruning result. How many SNPs survived the LD pruning
process and how many were lost?

<br>

### Generate an LD-pruned SNP list

We will use R to generate an LD-pruned SNP list in a format that can be
used by ANGSD for downstream analyses. This can be run **on the AWS
server**.

To launch R on the AWS server, run the following command.

``` bash
conda deactivate
conda activate lcwgs
R
```

Once you have launched R, run the following script in R to generate the
SNP list.

``` r
library(tidyverse)
basedir="~/day3/"
## position of SNPs that survived LD pruning
pruned_position <- read_lines(paste0(basedir, "ngsld/MME_ANGSD_PCA_subsampled_unlinked.id")) %>%
  str_remove("Mme_chr24_2558528-4558528:") %>%
  as.integer()
## filter SNPs in the original SNP list to only include those that survived LD pruning
pruned_snp_list <- read_tsv(paste0(basedir, "angsd/MME_ANGSD_PCA.mafs.gz")) %>%
  dplyr::select(1:4) %>%
  filter(position %in% pruned_position)
write_tsv(pruned_snp_list, paste0(basedir, "ngsld/LDpruned_snps.list"), col_names = F)
```

Note that this is a simple example with only one chromosome. If multiple
chromosomes are considered, both the chromosome name and position in the
SNP list need to match to the LD pruning output.

<br>

## Practical considerations

LD estimation and pruning tend to be computationally expensive because
for a chromosome with n SNPs, there are n<sup>2</sup> pairs of SNPs. To
make this computationally tractable, it is important to:

- consider downsampling your SNP list.
- impose stringent filters on minor allele frequencies since rare
  alleles are not helpful for LD estimation or for many of the
  downstream analyses. This can be done with the `--min_maf` flag
  (e.g. `--min_maf 0.05`).
- take advantage of the `--max_kb_dist` and `--max_snp_dist` flags for
  LD estimation and the `--max_dist` flag in LD pruning to limit the
  pairs of SNPs the program needs to consider.

Below are some specific recommendations with regard to two LD-related
applications.

For LD pruning, we recommend you to estimate LD with a relatively large
maximum distance (i.e. a distance at which you are very confident that
LD has decayed to the background level). Then, you can visualize the
patterns of LD decay by plotting LD estimates against physical distance
and fit a smooth curve to it. From this plot, you will be able to find
the smallest distance at which LD decays to the background level. Use
this value as your `--max_dist` in LD pruning. If everything still takes
too long, you can consider running each chromosome in parallel if you
have access to a computer cluster, or downsample the SNP list as we did
for this exercise. Also, if there are putative inversions that are
segregating in your samples, you can just remove all the SNPs that fall
into the inversions from the SNP list before LD pruning, because LD
extends too far when inversions are present.

For the visualization of large LD blocks, we recommend you to
drastically downsample the SNP list as we did here, and impose a very
stringent minor allele frequency filter. You will need to use a maximum
distance value that is large enough to span the entire LD block, so
having too many SNPs will significantly slow things down. If you are
just zooming into one of the edges of LD blocks (e.g. to identify the
breakpoint), you may use a denser SNP list as the maximum distance value
would be much smaller.

<br>
