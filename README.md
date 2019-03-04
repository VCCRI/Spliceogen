# Spliceogen
Spliceogen is an integrative, scalable tool for the discovery of splice-altering variants. Variants are assessed for their potential to create or disrupt any of the cis motifs which guide splice site definition: donors, acceptors, branchpoints, enhancers and silencers. Spliceogen integrates predictions from MaxEntScan<sup>1</sup>, GeneSplicer<sup>2</sup>, ESRseq<sup>3</sup> and Branchpointer<sup>4</sup>. Spliceogen accepts standard VCF/BED inputs and handles both SNPs and indels.

Author and maintainer: Steve Monger - s.monger@victorchang.edu.au

## Getting Started
Instructions for installation and obtaining the required genome annotation files.
### Installation:
Navigate to your desired installation directory and clone this repository:
```
git clone https://github.com/VCCRI/Spliceogen.git Spliceogen
```
### Dependencies:
-Bedtools

### Required annotation files:
-Any whole genome FASTA (.fa)

-Any GTF genome annotation (.gtf)
### Downloading required files:
Browse and download FASTA/GTF versions from [Gencode](https://www.gencodegenes.org/human/)

Alternatively, some recent (as of 2018) hg38 releases can be retrieved using:
```
> wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
> gunzip Homo_sapiens.GRCh38.dna.alt.fa.gz
> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz
> gunzip gencode.v29.basic.annotation.gtf.gz
```
### Including Branchpointer:
To include Branchpointer predictions for SNPs, the package must first be installed from an R command prompt.
```
> source("https://bioconductor.org/biocLite.R")
> biocLite("branchpointer")
```
Currently, only the development version of Branchpointer supports indels. To install this version instead:
```
> library(devtools)
> install_github("betsig/branchpointer_dev")
```
Then to include Branchpointer predictions, include the flag -branchpointer

Or for branchpointer_dev which handles both SNPs and indels, include the flag -branchpointerIndels 

## Running Spliceogen

### Basic Usage:
```
> cd path/to/Spliceogen
> ./RUN.sh -inputVCF path/to/singleOrMultipleFiles.vcf -fasta path/to/hgXX.fa -gtf path/to/annotation.gtf
```
Small VCF, BED, GTF and FASTA files are provided to demonstrate input and output format. Run this small example using the following command:

```
> ./RUN.sh -inputVCF toy/toy.vcf -gtf toy/toy.gtf -fasta toy/toy.fa
```

### BED input:
For BED inputs, replace the -inputVCF flag with -inputBED. See toy.bed for an example input format.
## Output

### Column labels:

The following abbreviations are used in the output header:

mes = MaxEntScan

gs = GeneSplicer

don = Donor

acc = Acceptor

ref = Reference allele

alt = Alternative allele

ESS = Silencer (ESRseq)

ESE = Enhancer (ESRseq)

So for example, the column "gsDonRef" contains GeneSplicer scores representing donor motif strength for the reference sequence, whereas "mesDonAlt" consists of MaxEntScan scores representing acceptor motif strength for the alternative sequence.

### Output Files

All scores and predictions can be found in the Spliceogen/output directory in a tab delimited format suitable for ANNOVAR annotation. Multiple output files are provided for each input VCF/BED. This includes one master file containing all scores for all variants, as well as several additional files containing only variants identified as most likely to be disruptive, ranked in descending order. The specific files generated are as follows:

1) "$file"_out.txt:

    Contains all scores generated for every variant, sorted in ascending chromosomal/start position order.

2) "$file"_withinSS.txt:

    Contains all variants that overlap annotated splice sites, alongside relevant scores and gene/exon information. Variants are sorted by donor/acceptor score decrease, such that the variants most likely to disrupt existing donor/acceptor splice sites appear at the top of this file.

3) "$file"_donorCreating.txt

    Contains variants outside of existing splice sites that are predicted to create donor motifs, ranked by P value, based on a logistic regression model of the MaxEntScan and GeneSplicer scores of known donor creating variants (discussed below)

4) "$file"_acceptorCreating.txt

    Same as above, but for acceptor creating variants.

5) "$file"_bpOutput.txt

    Contains Branchpointer prediction scores, including whether the variant is predicted to create or remove a branchpoint, based on the recommended Branchpointer thresholds.

## Database
Genome-wide SNV databases are available for [download](https://github.com/VCCRI/Spliceogen/tree/master/database). These cover all possible SNVs at every genomic position within all protein-coding multi-exon transcripts (based on Gencode v29). Both hg19 and hg38 are available.

An extensive version of the database containing MaxEntScan, GeneSplicer and ESRseq prediction scores for every possible SNV is available.

Alternatively, for users only interested in identifying variants which disrupt or create donor/acceptor motifs (and not silencer/enhancer predictions), a much leaner version of the database is provided. This version includes donor/acceptor predictions for all SNVs within splice sites, as well as all SNVs outside of splice sites that are likely to create donor/acceptor motifs. The database size is vastly reduced by excluding the large majority of SNVs which do not overlap splice sites and are not predicted to create a donor/acceptor motif.

## References:
1. Yeo, G., Burge, C., "Maximum Entropy Modeling of Short Sequence Motifs with Applications to RNA Splicing Signals", J Comput Biol. 2004; 11(2-3):377-94

2. Pertea, M., Lin, X., Salzberg, S., "GeneSplicer: a new computational method for splice site prediction", Nucleic Acids Res. 2001; 29(5):1185-90

3. Shendong, K., et al., "Quantitative evaluation of all hexamers as exonic splicing elements", Genome Res. 2011; 21(8): 1360-1374

4. Signal, B., et al., "Machine learning annotation of human branchpoints", Bioinformatics. 2018; 34(6):920-927
