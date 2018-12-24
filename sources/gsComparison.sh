cd /g/data/a32/quarterly_x1_10TB/WGS/splicing/combined
inputFile=$1
echo "inputFile: $inputFile"
echo "running multi line test"
time bin/linux/genesplicerAdapted $inputFile human > temp/gsOutputBenchmark.txt
echo "running single line test"
time ./gsSingleAdapted.sh $inputFile
