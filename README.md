# Spliceogen
Spliceogen is an integrative, scalable tool for the discovery of splice-disrupting variants. Variants are assessed for their potential to create or disrupt any of the cis motifs which guide splice site definition: donors, acceptors, branchpoints, enhancers and silencers. Spliceogen integrates some of the individually best performing models for splice motif prediction: MaxEntScan<sup>1</sup>, GeneSplicer<sup>2</sup>, ESRseq<sup>3</sup> and Branchpointer<sup>4</sup>. Spliceogen accepts standard VCF/BED inputs and handles both SNPs and indels.

## Getting Started

Instructions for installation and obtaining the required genome annotation files.

### Installation:

Navigate to your desired installation directory and clone this repository:

```
git clone https://github.com/VCCRI/Spliceogen.git Spliceogen
```

### Required annotation files:
-Any whole genome fasta (.fa)

-Any GTF genome annotation (.gtf)

### Downloading required files:

Browse and download desired versions from [UCSC] (hgdownload.soe.ucsc.edu/downloads.html#human)
and [Gencode] (https://www.gencodegenes.org/human)

Alternatively, hg38 releases can be retrieved using:
```
> wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
> gunzip Homo_sapiens.GRCh38.dna.alt.fa.gz
> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz
> gunzip gencode.v29.basic.annotation.gtf.gz
```
The "basic" gencode datasets are recommended. Users may want to restrict the analysis to annotated protein coding genes. Modify the GTF annotation using:
```
grep "protein_coding"
```

## Basic Usage:

```
> cd path/to/Spliceogen
> ./RUN.sh -inputVCF path/to/singleOrMultipleFiles.vcf -fasta path/to/hgXX.fa -gtf path/to/annotation.gtf
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

### BED input:

For BED inputs, replace the -inputVCF flag with -inputBED

## References:

1. Yeo, G., Burge, C., "Maximum Entropy Modeling of Short Sequence Motifs with Applications to RNA Splicing Signals", J Comput Biol. 2004; 11(2-3):377-94

2. Pertea, M., Lin, X., Salzberg, S., "GeneSplicer: a new computational method for splice site prediction", Nucleic Acids Res. 2001; 29(5):1185-90

3. Shendong, K., et al., "Quantitative evaluation of all hexamers as exonic splicing elements", Genome Res. 2011; 21(8): 1360-1374

4. Signal, B., et al., "Machine learning annotation of human branchpoints", Bioinformatics. 2018; 34(6):920-927
