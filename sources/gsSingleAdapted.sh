#!/bin/bash
inputFile=$(echo $1)
while read line1 && read line2; do
    echo -e "$line1\n$line2" > "$inputFile"SingleLine.FASTA
    output=$(bin/linux/genesplicerAdapted "$inputFile"SingleLine.FASTA human)
    if [ ! -z "$output" ]; then
        echo "$output" >> "$inputFile"gsScoresSingle.txt
    fi
done <"$inputFile"
