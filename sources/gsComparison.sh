#!/bin/bash
inputFile=$1
echo "inputFile: $inputFile"
echo "running multi line test"
time bin/linux/genesplicerAdapted $inputFile human > temp/gsOutputBenchmark.txt
echo "running single line test"
singleLine() {
    while read line1 && read line2; do
        echo -e "$line1\n$line2" > "$inputFile"SingleLine.FASTA
        output=$(bin/linux/genesplicerAdapted "$inputFile"SingleLine.FASTA human)
        if [ ! -z "$output" ]; then
            echo "$output" >> "$inputFile"gsScoresSingle.txt
        fi
    done <"$inputFile"
}
time singleLine
