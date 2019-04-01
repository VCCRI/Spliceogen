We provide two versions of the Spliceogen database. Both databases have genome-wide coverage, assessing every SNV at every position within every annotated multi-exon protein-coding transcript (1.29 billion base pairs in total, or 4.9 billion SNVs). They are available for both hg19 and hg38.

The “focussed” version contains all donor and acceptor predictions: 

hg19- https://s3-us-west-2.amazonaws.com/spliceogen/databases/hg19_focussed.zip

hg38- https://s3-us-west-2.amazonaws.com/spliceogen/databases/hg38_focussed.zip

The comprehensive version contains all donor, acceptor, silencer and enhancer predictions:

hg19- https://s3-us-west-2.amazonaws.com/spliceogen/databases/hg19.zip

hg38- https://s3-us-west-2.amazonaws.com/spliceogen/databases/hg38.zip

The focussed database contains predictions for all SNVs within annotated splice sites and all SNVs that are likely to create a de novo donor or acceptor motif. By excluding the vast majority of SNVs which fall outside of splice sites and are unlikely to create a donor/acceptor motif (logistic regression prediction score <0.7), this database is massively reduced in size without reducing the sensitivity of its donor/acceptor predictions.

Due to the sheer number of scores and predictions provided, we expect that the comprehensive database may be unwieldy for many use cases. In general we recommend running the tool to obtain comprehensive predictions, which has the advantage of including predictions for indels and (optionally) branchpoints, and the flexibility of selecting/customising your GTF annotation.
