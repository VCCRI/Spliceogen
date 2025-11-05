# Spliceogen

### Installation:
Clone this branch:
```
git clone -b comp3900 https://github.com/VCCRI/Spliceogen.git
```

### Spliceogen Dependencies:
Run ./setup.py to install Java and Bedtools. I just tested this on a fresh Ubuntu 24.04 VM. Everything installed and Spliceogen ran correctly.

## Running Spliceogen

First, try to run the "test usage" command. If you get the following JSON output, then the software is running correctly:


{"#CHR":"chr1", "START":"65510", "END":"65510", "REF":"T", "ALT":"G", "GENE":".", "withinSite":".", "mesDonRef":".", "mesDonAlt":".", "mesAccRef":"9.11", "mesAccAlt":"7.79", "gsDonRef":".", "gsDonAlt":".", "gsAccRef":"4.633193", "gsAccAlt":"2.350955", "ESEmaxRef":".", "ESEmaxAlt":".", "ESSminRef":"-0.286", "ESSminAlt":"-0.355", "donGainP":"0.02", "accGainP":"0.74", "donLossP":".", "accLossP":"."}


Next, setup the GTF and FASTA files, which will allow running any variant in the genome.

### Test Usage:
```
> cd path/to/Spliceogen
> ./RUN_helper.sh -input chr1:65510:T:G -gtf toy/toy.gtf -fasta toy/toy.fa
```

### Required whole genome files:

```
> wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
> gunzip GRCh38.primary_assembly.genome.fa.gz
> wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
> gunzip gencode.v49.basic.annotation.gtf.gz
```

### Usage:
```
> cd path/to/Spliceogen
> ./RUN_helper.sh -input chr1:65510:T:G -gtf <path>/gencode.v49.basic.annotation.gtf -fasta <path>/GRCh38.primary_assembly.genome.fa
```
