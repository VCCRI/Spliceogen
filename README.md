# Spliceogen

### Installation:
Clone this branch:
```
git clone -b comp3900 https://github.com/VCCRI/Spliceogen.git
```

### Spliceogen Dependencies:
-Bedtools

-Java

See setup.py. I ran those commands on a fresh Ubuntu 22.04 VM. Everything installed and Spliceogen ran correctly.

## Running Spliceogen

First, try to run the "test usage" variant. To run any variant in the genome, you will need the large annotation files below.

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
