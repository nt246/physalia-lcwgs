#### 3. Genotype calling

Now that you know how to calculate allele frequencies, we'll use these estimates to obtain prior probabilities 
for genotypes, which ultimately allows us to calculate genotype posterior probabilities and call genotypes.

In ANGSD, the option to call genotypes is `-doGeno`:

```bash
angsd -doGeno
...
-doGeno 0
        1: write major and minor
        2: write the called genotype encoded as -1,0,1,2, -1=not called
        4: write the called genotype directly: eg AA,AC etc
        8: write the posterior probability of all possible genotypes
        16: write the posterior probability of called genotype
        32: write the posterior probabilities of the 3 gentypes as binary
        -> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3
        -postCutoff=0.333333 (Only genotype to missing if below this threshold)
        -geno_minDepth=-1       (-1 indicates no cutof)
        -geno_maxDepth=-1       (-1 indicates no cutof)
        -geno_minMM=-1.000000   (minimum fraction af major-minor bases)
        -minInd=0       (only keep sites if you call genotypes from this number of individuals)

        NB When writing the posterior the -postCutoff is not used
        NB geno_minDepth requires -doCounts
        NB geno_maxDepth requires -doCounts
```

Therefore, if we set `-doGeno 2`, genotypes are coded as '0', '1', or '2' minor alleles. If we wanted to instead print the posterior 
probability of the genotype we could use '-doGeno 16'. The numeric arguments to `-doGeno` are additive so if we wanted to print *both*
the called genotype as 0, 1, or 2 *and* the probability of that genotype, we could use `-doGeno 18` (2 + 16). If we also wanted to 
print the major and minor allele we would use `-doGeno 19` (2 + 16 + 1).

To calculate the posterior probability of genotypes we need to specify a prior probability.
```
angsd -doPost
...
-doPost 0       (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
        4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an
...

```

We'll calculate genotype posterior probabilities using a HWE prior (`-doPost 1`) based on the allele frequencies estimated with `-doMaf 1`
and then output the posterior probabilities for the major/major, major/minor, minor/minor genotypes for each individual with `-doGeno 8` for
the PANY population. We'll limit our analysis to PANY biallelic SNPs (`SNP_pval 1e-6`). Note that while we need to estimate allele frequencies we already 
have a file containing them, which we don't need to write again (I/O can be expensive), so we can suppress writing another maf file by making the value for 
`-doMaf` negative. We'll use the BAM files as input (meaning that we have to recalculate GLs with `-GL`). You could also use pre-calculated genotype likelihoods
as input with `-glf` (binary) or `-glf10_text` (text). We'll use many of the same quality controls that we've been using throughout.

Note that while we need to estimate allele frequencies we already have a file containing them, which we don't need to write again (I/O can be expensive). In more recent
versions of ANGSD you can suppress writing another maf file by making the value for `-doMaf` negative, e.g. `-doMaf -1`. But for the version we have to write 
the maf file again.

```
angsd -b $DIR/PANY_bams.txt -ref $REF -out $RESDIR/PANY \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
   -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -rmTriallelic 0.05 -doPost 1 -doGeno 8
```
The first two columns of the output are the chromosome and position. The following columns list the the posterior probabilites
for the major/major, major/minor, and minor/minor for each individual in the same order as they appeared from the top of the bam list.

Take a look at the ouput:
```bash
less -S $RESDIR/PANY.geno.gz
```

**QUESTION**

What are the three posterior probabilites for PANY_07 at Mme_chr24:2558528-4558528 459780? Think back to how you extracted information 
from the glf files.

<details>

<summary> click for help </summary>

```bash
# find position of PANY_07 in the BAM file

INDNUM=$(grep -n "PANY_07.bam$" $DIR/PANY_bams.txt | cut -f1 -d':')
echo "$INDNUM"
```

So PANY_07 is at row 7 of the BAM list. Now extract their genotype probabilities from the .geno file for Mme_chr24:2558528-4558528 459780.

```bash
zcat $RESDIR/PANY.geno.gz | grep -m 1 $'^Mme_chr24:2558528-4558528\t459780\t' | cut -f 3- | perl -se '$start=($n-1)*3; @arr = split(/\t/,<>); print "@arr[$start .. $start+2]\n"' -- -n=$INDNUM
```

The most probably genotype configuration is major/major with a posterior probability of 0.984109.

</details> 

What do the genotype probabilities for PANY_07 at Mme_chr24:2558528-4558528 459780 sum to?

<details>

<summary> click for answer </summary>

1

</details>

Repeat this for PANY_03. What are their genotype probabilites are site Mme_chr24:2558528-4558528 459780?

<details>

<summary> Click for help </summary>

```bash
# Extract row for PANY_03
INDNUM=$(grep -n "PANY_03.bam$" $DIR/PANY_bams.txt | cut -f1 -d':')

# Extract the PANY_03's genotype posterior probabilities
zcat $RESDIR/PANY.geno.gz | grep -m 1 $'^Mme_chr24:2558528-4558528\t459780\t' | cut -f 3- | perl -se '$start=($n-1)*3; @arr = split(/\t/,<>); print "@arr[$start .. $start+2]\n"' -- -n=$INDNUM
```
The genotype posterior probabilities are 0.631517 0.326327 0.042156.

</details>

**EXERCISE**

Calculate the posterior probabilities for the PANY samples again, but this time using a uniform prior.

<details>

<summary> Click for help </summary>

```bash
angsd -b $DIR/PANY_bams.txt -ref $REF -out $RESDIR/PANY_unif \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
   -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -rmTriallelic 0.05 -doPost 2 -doGeno 8
```

</details>


Now extract the genotype probabilities for PANY_07 Mme_chr24:2558528-4558528 459780 from this new .geno file. What are the 
genotype probabilities?

<details>

<summary> Click for help </summary>

```bash
# find position of PANY_07 in the BAM file
INDNUM=$(grep -n "PANY_07.bam$" $DIR/PANY_bams.txt | cut -f1 -d':')

# Extract the genotype probablities
zcat $RESDIR/PANY_unif.geno.gz | grep -m 1 $'^Mme_chr24:2558528-4558528\t459780\t' | cut -f 3- | perl -se '$start=($n-1)*3; @arr = split(/\t/,<>); print "@arr[$start .. $start+2]\n"' -- -n=$INDNUM
```
The three genotype posterior probabilities are 0.969698 0.030302 0.000000.

</details>

Now extract the genotype posterior probabilities for PANY_03 at me_chr24:2558528-4558528 459780.

<details>

<summary> Click for help </summary>

```bash
# find position of PANY_03 in the BAM file
INDNUM=$(grep -n "PANY_03.bam$" $DIR/PANY_bams.txt | cut -f1 -d':')

# Extract the genotype probablities
zcat $RESDIR/PANY_unif.geno.gz | grep -m 1 $'^Mme_chr24:2558528-4558528\t459780\t' | cut -f 3- | perl -se '$start=($n-1)*3; @arr = split(/\t/,<>); print "@arr[$start .. $start+2]\n"' -- -n=$INDNUM
```

The genotype posterior probabilites are 0.333333 0.333333 0.333333.

</details>


**QUESTION**

What do the posterior genotype probabilities of `0.333333 0.333333 0.333333` under the uniform prior mean?

<details>

<summary> Click for answer </summary>

It means that the individual has missing data. There is no way of knowing what their genotype is or could be 
under an uninformative prior so the same probability is assigned to every genotype.

</details>

**QUESTION**

What happens when individuals having missing data when using a HWE prior `-doPost 1`?

<details>

<summary> Click for answer </summary>

Each possible genotype is assigned it's probability under the Hardy-Weinberg model based on the minor 
allele frequency estimated with `-doMaf`.

</details>


Now we'll call genotypes based on their maximum posterior probability and output the genotype in 0, 1, 2 format. We 
can set the genotype for an indivdiual to missing (-1) if the maximum posterior probability is less than a certain value. 
In this example we'll use a cutoff of 0.95. We'll also go ahead and write the major and minor alleles as well.

```bash
angsd -b $DIR/PANY_bams.txt -ref $REF -out $RESDIR/PANY_call \
   -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
   -minMapQ 20 -minQ 20 -minInd 5 -setMinDepthInd 1 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
   -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -rmTriallelic 0.05 -doPost 1 -doGeno 3 -postCutoff 0.95
```
Have a look at the file:
```bash
less -S $RESDIR/PANY_call.geno.gz
```

How many sites are in the PANY_call.geno.gz file?

```bash
zcat $RESDIR/PANY_call.geno.gz  | wc -l
```

How many sites have at least one individual with missing data?
```bash
zcat $RESDIR/PANY_call.geno.gz | grep -c "\-1"
```

How many sites have at least 10 individuals with missing data?

```bash
zcat $RESDIR/PANY_call.geno.gz | perl -ne '$n = () = $_ =~ /\t\-1/g; print "$n\n"' | grep -c "10$"
```

**Exercise**

Write a for loop to to calculate how many sites have missing data for 1,2,3,...,15 individuals.

<details>

<summary> Click for help </summary>

```bash
for NMISSING in {1..15};
do
	if [ $NMISSING == 1 ]; then printf "%s\t%s\n" 'NUMBER_MISSING_IND' 'NUMBER_SITES'; fi
	printf "%d\t%d\n" $NMISSING `zcat $RESDIR/PANY_call.geno.gz | perl -ne '$n = () = $_ =~ /\t\-1/g; print "$n\n"' | grep -c "$NMISSING$"`
done
```

</details>

Now try calling genotypes using a minimum posterior probability cutoff of 0.90 (and maybe some other cutoffs - have fun!). 
How many sites are in the file now? 

Calculate a spectrum for the number of individuals with uncalled genotypes as above. What do you make of this?

In general the `postCutoff` threshold will depend on the mean sequencing depth of your data as well as your application, i.e. 
how much miscalled genotypes could throw off an analysis.

Now you know how to
* Impose data quality filtering parameters
* Calculate genotype likelihoods
* Estimate allele frequencies and call SNPs based on these frequencies
* Incorporate allele frequencies into a genotype prior to calculate genotype posterior probabilities
* Call genotypes based on posterior probabilities

This concludes the Day 2 exercises.