LD estimation
================

  - [Data preparation](#data-preparation)
  - [Estimate LD](#estimate-ld)
  - [Visualize LD estimation result](#visualize-ld-estimation-result)
  - [Estimate LD decay](#estimate-ld-decay)
  - [LD pruning](#ld-pruning)

## Data preparation

## Estimate LD

``` bash
BASEDIR=/workdir/physalia-lcwgs/day_3/
NGSLD=/workdir/programs/ngsLD/
## Prepare a position file
zcat $BASEDIR/angsd/MME_ANGSD_PCA.mafs.gz | cut -f 1,2 | tail -n +2 | gzip > $BASEDIR/ngsld/MME_ANGSD_PCA.pos.gz
## Run ngsLD
$NGSLD/ngsLD \
--geno $BASEDIR/angsd/MME_ANGSD_PCA.beagle.gz \
--pos $BASEDIR/ngsld/MME_ANGSD_PCA.pos.gz \
--n_ind 60 \
--n_sites 865 \
--out $BASEDIR/ngsld/MME_ANGSD_PCA.ld \
--probs \
--rnd_sample 1 \
--seed 1 \
--n_threads 10
```

## Visualize LD estimation result

``` bash
BASEDIR=/workdir/physalia-lcwgs/day_3/
NGSLD=/workdir/programs/ngsLD/
cat $BASEDIR/ngsld/MME_ANGSD_PCA.ld | bash $NGSLD/scripts/LD_blocks.sh Mme_chr24 2558408 2564908
```

![](../ngsld/LD_blocks.r2.pdf)

## Estimate LD decay

``` bash
Rscript --vanilla --slave $NGSLD/scripts/fit_LDdecay.R --ld_files ld_files.list --out plot.pdf
```

``` r
library(tidyverse)
ld <- read_tsv("../ngsld/MME_ANGSD_PCA.ld", col_names = F)
ld %>% 
  ggplot(aes(x=X3, y=X7)) +
  geom_point() +
  geom_smooth()
```

## LD pruning

``` bash
BASEDIR=/workdir/physalia-lcwgs/day_3/
NGSLD=/workdir/programs/ngsLD/
perl $NGSLD/scripts/prune_graph.pl --in_file $BASEDIR/ngsld/MME_ANGSD_PCA.ld --max_kb_dist 5 --min_weight 0.5 --out $BASEDIR/ngsld/MME_ANGSD_PCA_unlinked.id
```
