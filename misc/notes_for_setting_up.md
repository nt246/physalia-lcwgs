Notes for setting up
================

## AWS log in

```bash
chmod 600 Downloads/lcwgs24.pem
ssh -i Downloads/lcwgs24.pem ubuntu@34.218.244.188
```

## Software installation

#### lcwgs conda environment

* Most software packages for this tutorial can be installed into a single conda environment without experiencing any conflicts.  
* These include fastqc, trimmomatic, picard, samtools, bowtie2, bamutil, r, tidyverse  
* Carlo installed most of these in 2024 into a conda environment called `lcwgs`. I have generated [a yaml file](lcwgs.yaml) from this environment so that it can be created again. 

#### GATK

* GATK is a bit more complicated, especially since we need (an older version (v3.7))[https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2;tab=live_object]. 
* It may be easier installed it in a separate conda environment to avoid conflicts.
* It is also required to download the source file from the Internet, upload it to the server, unzip it, and provide its path to the gatk command.  

```bash
## create a separate conda environment
mamba create -c bioconda -c conda-forge -n gatk-3.7 gatk=3.7
## download GATK-3.7 from https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2;tab=live_object
## upload from local to AWS (paths need to be changed accordingly)
scp -i Downloads/lcwgs24.pem  Downloads/package-archive_gatk_GenomeAnalysisTK-3.7-0-gcfedb67.tar ubuntu@34.217.81.210:~/Share/gatk-3.7/
## unzip on AWS
cd ~/Share/gatk-3.7/
tar -xvf package-archive_gatk_GenomeAnalysisTK-3.7-0-gcfedb67.tar
## provide the path the the .jar file to gatk
conda activate gatk-3.7
gatk-register ~/Share/gatk-3.7/GenomeAnalysisTK.jar
```

The `gatk` command is then available in the `gatk-3.7` conda environment

#### ngsLD

* ngsLD is also rather complicated because it isn't available on conda and has a lot of dependencies.
* It may be easier to create a separate conda environent with all required dependencies, and install ngsLD manually from there
* Note that LDHeatmap can be installed to this conda environment, so that students don't have to install it themselves (which tends to be tricky).

```bash
## create a separate conda environment with required dependencies
mamba create -c bioconda -c conda-forge -c genomedk -n ngsLD gcc zlib gsl pandas graph-tool r-optparse r-ggplot2 r-reshape2 r-plyr r-gtools r-ldheatmap
## the installation step had problems finding these following paths, so I had to specify them. Note that the last two lines need to be run first before running ngsLD after each login.
export CPATH=/home/ubuntu/src/conda/envs/ngsLD/include/
export PKG_CONFIG_PATH=/home/ubuntu/src/conda/envs/ngsLD/lib/pkgconfig/:$LD_LIBRARY_PATH
export LIBRARY_PATH=/home/ubuntu/src/conda/envs/ngsLD/lib/
export LD_LIBRARY_PATH=/home/ubuntu/src/conda/envs/ngsLD/lib/:$LD_LIBRARY_PATH
## download and install ngsLD from GitHub and compile
cd ~/Share
git clone https://github.com/fgvieira/ngsLD.git
make
```

Again, note that the `LIBRARY_PATH` and `LD_LIBRARY_PATH` variables need to be specified after each login.

#### NGSAdmix

NGSAdmix is included in ANGSD, so it can be invoked from wherever ANGSD installed, or by having ANGSD in a conda environment. We took the former approach this year, with `NGSADMIX=/home/ubuntu/angsd/misc/NGSadmix`

#### PCAngsd

Use the following code to install PCANngsd in a conda environment.

```bash
cd Share/
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
mamba env create -f environment.yml
conda activate pcangsd
python setup.py build_ext --inplace
pip3 install -e . ## if you run into issues pertaining to ssl certificates, try "pip3 install --trusted-host pypi.org -e ." instead
conda deactivate
```

## Directory setup

#### Day 1

*Use the following code to prepare the directory structure for students to copy and paste

```bash
cd ~/Share
git clone https://github.com/nt246/physalia-lcwgs.git
cp -r physalia-lcwgs/day_1/ day1/
rm -rf day1/img/
rm -rf day1/markdowns/
rm -f day1/fastqc/*
mkdir day1/adapter_clipped
mkdir day1/bam
```

#### Day 2

* Note that the renamed refrence genome and bam files do not yet exist. It would be helpful to have these files on GitHub.

#### Day 3

* Students can get the files that they need directly from the GitHub repo by following to the tutorial. No preparation needed. 