#!/bin/bash
#assuming hg38 for now
grep "^#" toy/toy.vcf > inputString.vcf
echo "$2" | awk -F'[:]' -v OFS="\t" '{print $1, $2, ".", $3, $4, ".", ".", "."}' >> inputString.vcf

rm output/inputString.vcf* 2>/dev/null

./RUN.sh "$@" > /dev/null 2>&1

#JSON output
if [ -f "output/inputString.vcf_out.txt" ]; then
    awk 'BEGIN{FS="\t"}
    NR==1 {for(i=1;i<=NF;i++) header[i]=$i; next}
    {
        printf "{";
        for(i=1;i<=NF;i++){
            printf "\"%s\":\"%s\"", header[i], $i;
            if(i<NF) printf ", "
        }
        print "}"
    }' output/inputString.vcf_out.txt
fi
