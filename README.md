# Spliceogen
Spliceogen is an integrative, scalable tool for the discovery of splice-altering variants. Variants are assessed for their potential to create or disrupt any of the cis motifs which guide splice site definition: donors, acceptors, branchpoints, enhancers and silencers. Spliceogen integrates scores from MaxEntScan<sup>1</sup>, GeneSplicer<sup>2</sup>, ESRseq<sup>3</sup> and Branchpointer<sup>4</sup>, and provides predictions based on logistic regression models trained on reported splice-altering variants<sup>5</sup>. Spliceogen accepts VCF/TSV inputs and handles both SNVs and indels.

Publication: https://doi.org/10.1093/bioinformatics/btz263

Maintainer: Steve Monger - s.monger@victorchang.edu.au

## Getting Started

### Installation:
Navigate to your desired installation directory and clone this repository:
```
git clone https://github.com/VCCRI/Spliceogen.git Spliceogen
```
### Required annotation files:
-Any whole genome FASTA (.fa)

-Any GTF genome annotation (.gtf)

### Obtaining required files:
FASTA/GTF files can be downloaded from [Gencode](https://www.gencodegenes.org/human/)

Alternatively, some recent (as of 2019) hg38 releases can be retrieved using:
```
> wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz
```
### Spliceogen Dependencies:
-Bedtools

-Java

## Running Spliceogen

### Basic Usage:
```
> cd path/to/Spliceogen
> ./RUN.sh -input path/to/variant/file(s) -fasta path/to/hgXX.fa -gtf path/to/annotation.gtf
```
Example VCF, TSV, GTF and FASTA files are provided to demonstrate input and output formats. Run this small example using the following command:

```
> ./RUN.sh -input toy/toy.vcf -gtf toy/toy.gtf -fasta toy/toy.fa
```

To include optional Branchpointer predictions, see instructions [below](#Including%20Branchpointer%3A)

### Input formats:

As an alternative to VCF, a custom tab-separated format is allowed (chr    start    ref    alt). Gzipped GTF/VCF/TSV files are accepted.

## Scalability

Spliceogen is highly scalable. Predictions are generated at a rate of 2.3 million variants/compute hour, with peak memory usage less than 500Mb. Benchmarking was performed using a single compute node with 1 CPU allocated, without Branchpointer predictions. If preferred, a genome-wide database of pre-computed predictions for all SNVs within genes is available. Contact us to obtain this.

## Output

### Predictions

We developed logistic regression models for each of the following classes of splice-altering variants: donor loss, acceptor loss, donor gain and acceptor gain. Using these, we derive probability values which are used to rank variants based on the likelihood that they will cause each kind of splice-altering variant. Variants outside of splice sites are assigned donor and acceptor gain scores only, while variants within donor/acceptor splice sites are assigned only donor/acceptor loss scores.  


Note that these probability values are used for ranking, and should not be interpreted as the actual probability of splice alteration. Similarly, these scores should not be compared across different models (for instance, donor loss and donor gain).

### Files

Multiple output files are created for each VCF/TSV in the Spliceogen/output directory. A master "\_out" file contains all scores for all variants, in a format suitable for ANNOVAR<sup>6</sup> variant annotation. Several additional files show predictions for variants identified as most likely to be disruptive, ranked in descending order. The specific files generated are as follows:

1) "$file"_out.txt:

    Contains all scores generated for every variant, sorted in standard ascending chromosomal/start order.

2) "$file"_withinSS.txt:

    Contains all variants that overlap annotated splice sites. The overlapping splice sites are denoted by their exonID and "\_donor" or "\_acceptor". Variants are sorted by the maximum of donLossP and accLossP, such that variants most likely to disrupt acceptor/donor splice sites appear at the top of this file.

3) "$file"_ssGain.txt

    Contains variants outside of existing splice sites that are predicted to create donor or acceptor motifs. Variants are sorted by the maximum of donGainP and accGainP, such that variants most likely to create acceptor/donor splice sites appear at the top of this file.

4) "$file"_bpOutput.txt

    Contains Branchpointer prediction scores, including whether the variant is predicted to create or remove a branchpoint, based on the recommended Branchpointer thresholds.

### Column labels:

The following abbreviations are used in the output headers:

donGainP = donor creation logistic regression probability value

accGainP = acceptor creation logistic regression probability value

donLossP = donor disruption logistic regression probability value

accLossP = acceptor disruption logistic regression probability value

withinSS = within splice site

don = Donor

acc = Acceptor

ref = Reference allele

alt = Alternative allele

mes = MaxEntScan

gs = GeneSplicer

ESS = exonic splicing silencer (ESRseq score)

ESE = exonic splicing enhancer (ESRseq score)

So for example, the column "gsDonRef" contains GeneSplicer scores representing donor motif strength for the reference sequence, whereas "mesDonAlt" consists of MaxEntScan scores representing acceptor motif strength for the alternative sequence.

## Including Branchpointer:

To include Branchpointer predictions, include the -branchpointer flag and specify the genome build:

```
*basic usage command* -branchpointer hgXX
```
Or for branchpointer_dev which handles both SNPs and indels, use the flag -branchpointerIndels hgXX

### Branchpointer dependencies:

-R (tested on v3.4.3)

-Branchpointer

-BSgenome

The current Bioconductor release of Branchpointer supports SNV predictions. To install it from an R prompt:

```
> source("https://bioconductor.org/biocLite.R")
> biocLite("branchpointer")
```
The development version of Branchpointer also supports indels. To install this version instead:
```
> library(devtools)
> install_github("betsig/branchpointer_dev")
```
From an R prompt, install the hg38 BSgenomes package using the below command. For hg19, edit the 2nd line to "hg19".

```
> source("https://bioconductor.org/biocLite.R")
> biocLite("BSgenome.Hsapiens.UCSC.hg38")
```

## References:
1. Yeo, G., Burge, C., "Maximum Entropy Modeling of Short Sequence Motifs with Applications to RNA Splicing Signals", J Comput Biol. 2004; 11(2-3):377-94

2. Pertea, M., Lin, X., Salzberg, S., "GeneSplicer: a new computational method for splice site prediction", Nucleic Acids Res. 2001; 29(5):1185-90

3. Shendong, K., et al., "Quantitative evaluation of all hexamers as exonic splicing elements", Genome Res. 2011; 21(8):1360-1374

4. Signal, B., et al., "Machine learning annotation of human branchpoints", Bioinformatics. 2018; 34(6):920-927

5. Shiraishi, Y., et al., "A comprehensive characterization of cis-acting splicing-associated variants in human cancer", Genome Res. 2018; 28(8):1111-1125
