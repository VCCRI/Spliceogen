#!/bin/bash
echo "$2" | tr ":" "\t" > inputString.tsv
rm output/inputString.tsv*

./RUN.sh "$@" > /dev/null 2>&1

#JSON output
if [ ! -f $FASTAPATH ]; then
    awk 'BEGIN{FS="\t"}
    NR==1 {for(i=1;i<=NF;i++) header[i]=$i; next}
    {
        printf "{";
        for(i=1;i<=NF;i++){
            printf "\"%s\":\"%s\"", header[i], $i;
            if(i<NF) printf ", "
        }
        print "}"
    }' output/inputString.tsv_out.txt
fi
