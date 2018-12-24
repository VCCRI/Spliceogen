#!/bin/bash
inputFile=$(echo $1)
numLines=$(cat "$inputFile" | wc -l)
numVariants=$(( $numLines / 2))
for i in $(seq 1 $numVariants); do
    #split into one fasta/file
    x=$(( i * 2))
    #fileID=$( head -${x} "$inputFile" | tail -2 | head -1 )
    head -${x} "$inputFile" | tail -2 > temp/gsInput"$fileID".FASTA
    output=$(bin/linux/genesplicerAdapted temp/gsInput"$fileID".FASTA human)
    if [ ! -z "$output" ]; then
        echo "$output" >> "$inputFile"gsScoresSingle.txt
    fi
done
