
## 2. Genotype likelihoods

Many of the calculations that you perform with ANGSD are based on **genotype likelihoods** (GL)s, which improves accuracy over hard-called genotypes when depth is low 
(see for instance [this paper](https://doi.org/10.1534/genetics.113.154740)). 
We will now calculate GLs for each individual at all sites in our dataset. 

To do so you need to specify which genotype likelihood model to use.
```
$angsd -GL
...
-GL=0: 
	1: SAMtools
	2: GATK
	3: SOAPsnp
	4: SYK
	5: phys
	6: Super simple sample an allele type GL. (1.0,0.5,0.0)
	7: outgroup gls
	-trim		0		(zero means no trimming)
	-tmpdir		angsd_tmpdir/	(used by SOAPsnp)
	-errors		(null)		(used by SYK)
	-minInd		0		(0 indicates no filtering)

Filedumping:
	-doGlf	0
	1: binary glf (10 log likes)	.glf.gz
	2: beagle likelihood file	.beagle.gz
	3: binary 3 times likelihood	.glf.gz
	4: text version (10 log likes)	.glf.gz
```
A brief description of these different models can be found [here](http://www.popgen.dk/angsd/index.php/Genotype_likelihoods).
The GATK model is based on the one from the first GATK paper, SAMtools uses a more sophisticated error model (does not assume errors are independent), 
SOAPsnp requires a reference sequence for recalibration of quality scores, SYK is error-type specific.
In most cases, GATK and SAMtools models should give similar results.

Let's first work with only the PANY samples.
A possible command to calculate genotype likelihoods is

```
$angsd -b $DIR/PANY_bams.txt -ref $REF -out $RESDIR/PANY \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 1 -doGlf 4
```

where we specify:
* -GL 1: genotype likelihood model as in SAMtools
* -doGlf 4: output in text format

Note that you can usually speed up computation with `-nThreads INT` when working with real data.

The `-C` option adjusts mapping quality for reads with excessive mismatches. The [samtools](http://www.htslib.org/doc/samtools-mpileup.html)
documentation recommends setting this to 50 for reads mapped with BWA. `-C` and `-baq` (not used here, but see below) require that you supply the reference genome with `-ref`.
`-baq` (base alignment quality) is another option to consider in order to adjust the base quality scores around INDELS.
<br>
`-minInd X` in conjunction with `-minIndDepth Y` requires at least *X* individuals to be covered by at least *Y* (default 1) reads to keep a site. 
`-setMinDepth` and `-setMaxDepth` set minimum and maximum total site depth (depth summed across all individuals) limits, respectively.


**QUESTION**
Looking at the standard output. What percentage of the sites provided to ANGSD were actually retained after applying the filters?

<details>

<summary> click here for help </summary>

```bash
-> Total number of sites analyzed: 1800717
-> Number of sites retained after filtering: 1070747 
```

59% of sites were retained.

</details>

ANGSD always dumps a log file with information on how it was run. Check it out:

```bash
less $RESDIR/PANY.arg
```

Have a look at the GLs with `less -S $RESDIR/PANY.glf.gz`. The first two columns refer to the reference sequence (chromososome and position). Then you have 10 likelihoods
for all possible genotypes in the order AA, AC, AG, AT, CC, CG, CT, GG, GT, TT. This set of 10 likelihoods is repeated sequentially starting from
the left of the file for each individual in the row order of individuals in the BAM file. The values are log likelihood ratios 
scaled to the most likely genotype, so the most likely genotype has a value of 0.

**QUESTION**
We are analyzing 15 individuals so we should have 152 fields in the glf file. You should confirm this and try to print the likelihoods
for the individual called PANY_04 at position Mme_chr24:2558528-4558528:34213 (each bam file is named by the individual, i.e. <idividual ID>.bam). What is 
their most likely genotype? If you need help you can click below.

<details>

<summary> click for help extracting GL info </summary>

```bash
# Count number of columns and subtract 2 (chromosome and position fields) to get the number of likelihood values

echo "$(($(zcat $RESDIR/PANY.glf.gz | head -n1 | wc -w)-2))"
```

You should see that indeed there are 150 likelihood values. Now figure out what line PANY_04 is at in the bam list.

```bash
INDNUM=$(grep -n "PANY_04.bam$" $DIR/PANY_bams.txt | cut -f1 -d':')
echo "$INDNUM"
```

So this individual is at row 4 in the bam list. Now we can extract their likelihoods.

```bash
zcat $RESDIR/PANY.glf.gz | grep -m 1 $'^Mme_chr24:2558528-4558528\t34213\t' | cut -f 3- | perl -se '$start=($n-1)*10; @arr = split(/\t/,<>); print "@arr[$start .. $start+9]\n"' -- -n=$INDNUM
```
For PANY_04 the 5th likelihood is zero, corresponding to the genotype 'CC'.

</details>

**QUESTION**
If you were to carry out what you just did for PANY_09 at site Mme_chr24:2558528-4558528 34213 you would see that there GLs are 
`0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000`. What does this mean?

<details>

<summary> click for answer </summary>

If all GLs are zero this means that there is no data for this individual.

</details>

**BONUS QUESTION**
Try to output genotype likelihood files in binary format (be sure to change the ouput name to avoid overwriting). Which option should you use? Can you open these files?
Look at the file sizes of text vs binary format. Which one is smaller?

<details>

<summary> click for help </summary>

Use `-doGlf 1` to output genotype likelihoods in binary. So the full command would be.

```
$angsd -b $DIR/PANY_bams.txt -ref $REF -out $RESDIR/PANY_binary \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 1 -doGlf 1
```

The binary .glf file is not human readable (it will just be a bunch of gibberish if you try to look at it). The compressed binary .glf is 
slightly smaller (66M) compared to the compressed text version (69M). While this size difference is small for this test dataset it can be 
appreciable for larger datasets. 

</details>

**BONUS QUESTION**
Try to change some filtering options and record the number of entries in the final output file (remember to change the output name to avoid overwriting).

You have learned how to calculate and read genotype likelihood files.
Next you are going to learn how to estimate allele frequencies with ANGSD.

[click here](https://github.com/nt246/physalia-lcwgs/blob/main/day_2/markdowns/03_allele_frequencies.md) to move to the next session.

--------------------------------
