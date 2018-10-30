#!/bin/bash
#set default parameters
POSITIONAL=()
INPUTVCF="FALSE"
INPUTBED="FALSE"
INPUTFILES=""
USEBP="FALSE"
USEBPINDELS="FALSE"
ANNOTATION=""
FASTAPATH=""
#parse command line args
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -inputVCF)
    INPUTVCF="TRUE"
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -inputBED)
    INPUTBED="TRUE"
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -fasta)
    FASTAPATH="$2"
    shift
    shift ;;
    -gtf)
    ANNOTATION="$2"
    shift
    shift ;;
    -branchpointer)
    USEBP="TRUE"
    USEBPINDELS="FALSE"
    shift ;;
    -branchpointerIndels)
    USEBPINDELS="TRUE"
    USEBP="FALSE"
    shift ;;
    *)
    POSITIONAL+=("$1")
    shift ;;
esac
done
set - "${POSITIONAL[@]}"
#input error handling
if [ ! -f $FASTAPATH ]; then
    echo -e "Fasta file not found: use -fasta ./path/to/hgXX.fa\nExiting..."
    exit 1
elif [ ! -f "$ANNOTATION" ]; then
    echo "GTF annotation file not found: use -annotation path/to/gencodeXX.gtf\nExiting..."
    exit 1
fi
#prepare splice site intervals from annotation.gtf
if [ ! -f data/annotationIntervals.txt ] || [[ "$ANNOTATION" -nt data/annotationIntervals.txt ]] ; then
    echo "Preparing splice site annotation..."
    grep '[[:blank:]]gene[[:blank:]]\|[[:blank:]]exon[[:blank:]]' data/"$ANNOTATION" | java -cp bin getSpliceSiteIntervalsFromGTF > data/annotationIntervals.txt
fi
#for each file
for FILE in $INPUTFILES; do
    fileID=$(echo "$FILE" | xargs -n 1 basename)
    echo "Input file: $fileID"
    #remove temp files from any previous run
    rm temp/"$fileID"* 2> /dev/null
    #sort body of input file
    grep "^#" "$FILE" > temp/"$fileID"_sorted
    grep -v "^#" "$FILE" | sort -k1,1 -k2,2n >> temp/"$fileID"_sorted
    #bedtools intersect to get strand info
    echo "Retrieving strand info..."
    grep '[[:blank:]]gene[[:blank:]]' data/"$ANNOTATION" | sort -k1,1 -k4,4n | bedtools intersect -a temp/"$fileID"_sorted -b stdin -wa -wb -sorted  > temp/"$fileID"unstrandedInput.txt
    #generate flanking intervals.bed for bedtools getfasta
    if [ ! "$INPUTVCF" = "FALSE" ]; then
        grep '[[:blank:]]+[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "+", $4, $5}' | sort -u | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        grep '[[:blank:]]-[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "-", $4, $5}' | sort -u | java -cp bin getFastaIntervals >> temp/"$fileID"fastaIntervals.bed
    elif [ ! "$INPUTBED" = "FALSE" ]; then
        grep '[[:blank:]]+[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "+", $5, $6}' | sort -u | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        grep '[[:blank:]]-[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "-", $5, $6}' | sort -u | java -cp bin getFastaIntervals >> temp/"$fileID"fastaIntervals.bed
    fi
    echo "Retrieving flanking FASTA sequence..."
    bedtools getfasta -fi $FASTAPATH -bed temp/"$fileID"fastaIntervals.bed -name -s > temp/"$fileID"seqToScan.FASTA
    #seqScan: generates input strings for maxentscan and genesplicer as well as ESRseq scores
    echo "Scanning for motifs..."
    java -cp bin seqScan temp/"$fileID"seqToScan.FASTA -useESR $fileID 1>&2
    #run maxentscan
    echo "Running MaxEntScan..."
    perl score5.pl temp/"$fileID"mesDonorInput.txt | java -cp bin processScoresMES > temp/"$fileID"mesDonorScores.txt
    retVal=( ${PIPESTATUS[0]} )
    if [ $retVal -ne 0 ]; then
        echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
    exit $retVal
    fi
    perl score3.pl temp/"$fileID"mesAcceptorInput.txt | java -cp bin processScoresMES > temp/"$fileID"mesAcceptorScores.txt
    retVal=( ${PIPESTATUS[0]} )
    if [ $retVal -ne 0 ]; then
        echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
    exit $retVal
    fi
    #run genesplicer
    echo "Running GeneSplicer..."
    bin/linux/genesplicerAdapted temp/"$fileID"gsInput.FASTA human > temp/"$fileID"gsScores.txt
    #run branchpointer SNPs
    if [ "$USEBP" = "TRUE" ]; then
        echo "Running Branchpointer..."
        Rscript --slave --vanilla bin/bpProcessing.R temp/"$fileID"input.txt $(pwd) "$BPGTF" &> output/"$fileID"bpLog.txt
        #awk -v OFS=\\t '{print $2, $3, $4, $8, $9, $15, $16, $21, $22, $23, $24}' output/bpOutputSNPs.txt > output/bpSNPsSummarised.txt
    fi
    #run branchpointer indels
    if [ "$USEBPINDELS" = "TRUE" ] ; then
        echo "Running Branchpointer..."
        while read -r na chr start strand ref alt; do
            refLength=${#ref}
            altLength=${#ref}
            end=$((refLength+start))
            if [ $altLength -gt 1 ] || [ $refLength -gt 1 ]; then
                echo -e ".\t$chr\t$start\t$end\t$strand\t$ref\t$alt" >> temp/"$fileID"bpInputIndels.txt
            fi
        done < temp/"$fileID"input.txt
        Rscript --slave --vanilla bin/bpProcessingINDELS.R temp/"$fileID"input.txt temp/"$fileID"bpInputIndels.txt $(pwd) "$BPGTF" &> output/"$fileID"bpLog.txt
        #awk -v OFS=\\t '{print $2, $3, $4, $8, $9, $16, $17, $23, $24, $25, $26}' output/bpOutputIndels.txt > output/bpIndelsSummarised.txt
        #awk -v OFS=\\t '{print $2, $3, $4, $8, $9, $15, $16, $22, $23, $24, $25}' $SNPpath"bpOutput_SNPs.txt" > output/bpSNPsSummarised.txt
    fi
    #merge scores into one line
    echo "Processing scores..."
    cat temp/"$fileID"mesDonorScores.txt temp/"$fileID"mesAcceptorScores.txt temp/"$fileID"gsScores.txt temp/"$fileID"ESRoutput.txt data/annotationIntervals.txt | sort -k1,1 -V -k 2,2n -k 3 -k 4 -s | java -cp bin mergeOutput > output/"$fileID"_out.txt
    #clean up temp files
    rm temp/"$fileID"* 2> /dev/null
done
