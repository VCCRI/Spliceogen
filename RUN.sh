#!/bin/bash
function printVersion {
    printf "Spliceogen 1.0 1-March-2019\n"
}
function printHelp {
    cat <<-END
Usage:
------
3 required args:
1)  -inputVCF path/to/VCF/input(s).VCF
        OR
    -inputBED path/to/input(s).BED 
        Note: wildcard matching of multiple files is allowed
2)  -gtf path/to/annotation.GTF
3)  -fasta path/to/genome.fasta
optional arg:
4)  -branchpointer hgXX
        OR
    -branchpointerIndels hgXX
        Note: user must specify hg19 or hg38
END
}
#set default parameters
POSITIONAL=()
INPUTVCF="FALSE"
INPUTBED="FALSE"
INPUTFILES=""
USEBP=""
USEBPINDELS=""
ANNOTATION=""
FASTAPATH=""
GENOMEBUILD=""
#parse command line args
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -h|--help)
    printHelp
    exit 0
    ;;
    -v|--version)
    printVersion
    exit 0
    ;;
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
    GENOMEBUILD="$2"
    USEBP="TRUE"
    USEBPINDELS="FALSE"
    shift
    shift ;;
    -branchpointerIndels)
    GENOMEBUILD="$2"
    USEBP="TRUE"
    USEBPINDELS="TRUE"
    shift
    shift ;;
    *)
    POSITIONAL+=("$1")
    shift ;;
esac
done
set - "${POSITIONAL[@]}"
#check input files exist and are not gzipped
if [ ! -f $FASTAPATH ]; then
    echo -e "Fasta file not found: use -fasta ./path/to/hgXX.fa\nExiting..."
    exit 1
elif [ ! -f "$ANNOTATION" ]; then
    echo "GTF annotation file not found: use -gtf path/to/gencodeXX.gtf\nExiting..."
    exit 1
fi
checkGzip=$( file --mime-type "$FASTAPATH" "$ANNOTATION" "$INPUTFILES" | grep gzip)
if [ ! "$checkGzip" = "" ]; then
    echo "Error: Input FASTA and GTF files must be unzipped. Exiting..."
    exit 1
fi
#check input files for consistenty in "chr" nomenlature
gtfChr=$(head "$ANNOTATION" | tail -1 | awk '{print $1}' | grep chr)
fastaChr=$(head -1 "$FASTAPATH" | awk '{print $1}' | grep chr)
firstInputFile=$(echo "$INPUTFILES" | awk '{print $1}')
inputChr=$(tail -1 "$firstInputFile" | awk '{print $1}' | grep chr)
warningChr="Warning: it appears the provided gtf, fasta, and input files use inconsistent Chromosome nomenclature. Eg. \"chr1\" vs \"1\". This will likely cause issues. Please edit them for consistency"
if [ "$gtfChr" == "" ]; then
    if [ "$fastaChr" != "" ]; then
        echo "$warningChr"
    elif [ "$inputChr" != "" ]; then
        echo "$warningChr"
    fi
else
    if [ "$fastaChr" == "" ]; then
        echo "$warningChr"
    elif [ "$inputChr" == "" ]; then
        echo "$warningChr"
    fi
fi
#prepare splice site intervals from annotation.gtf
gtfBasename=$(basename $ANNOTATION)
if [ ! -f data/"$gtfBasename"_SpliceSiteIntervals.txt ] || [[ "$ANNOTATION" -nt data/"$gtfBasename"_SpliceSiteIntervals.txt ]] ; then
    echo "Preparing splice site annotation..."
    grep '[[:blank:]]gene[[:blank:]]\|[[:blank:]]transcript[[:blank:]]\|[[:blank:]]exon[[:blank:]]' "$ANNOTATION" | grep -v '^GL000' |
    java -cp bin getSpliceSiteIntervalsFromGTF > data/"$gtfBasename"_SpliceSiteIntervals.txt
fi
#for each input VCF/BED file
for FILE in $INPUTFILES; do
    fileID=$(echo "$FILE" | xargs -n 1 basename)
    #check current file exists 
    if [ ! -f "$FILE" ]; then
        echo "Error: variant input file not found: $FILE \n Exiting..."
        exit 1
    fi
    echo "Input file: $fileID"
    #remove temp files from any previous run
    rm temp/"$fileID"* 2> /dev/null
    #sort body of input file
    grep "^#" "$FILE" > temp/"$fileID"_sorted
    grep -v "^#" "$FILE" | sort -k1,1 -k2,2n >> temp/"$fileID"_sorted
    #check bedtools is installed
    bedtoolsLocation=$(which bedtools);
    if [ "$bedtoolsLocation" == "" ]; then
        printf -- 'Warning: Bedtools does not appear to be installed.\n';
        printf -- 'Get it here: https://bedtools.readthedocs.io/en/latest/content/installation.html\n';
    fi;
    #bedtools intersect to get strand info
    echo "Retrieving strand info..."
    grep '[[:blank:]]gene[[:blank:]]' "$ANNOTATION" | sort -k1,1 -k4,4n | grep -v '^GL000' | bedtools intersect -a temp/"$fileID"_sorted -b stdin -wa -wb -sorted  > temp/"$fileID"unstrandedInput.txt 
    if [ $? -ne 0 ]; then
        echo "Warning. Bedtools intersect returned non-zero exit status. Intersection failed between provided variant VCF/BED file and provided GTF. See above error message for more details"
    fi
    if [ ! -s temp/"$fileID"unstrandedInput.txt ]; then
        echo "Error: no variants were returned following bedtools intersect between input file \""$fileID"\" and gtf. \n Exiting..."
        exit 1
    fi
    #generate flanking intervals.bed for bedtools getfasta and branchpointer input
    if [ "$INPUTVCF" = "TRUE" ]; then
        grep '[[:blank:]]+[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "+", $4, $5}' | ( [[ "$USEBP" ]] && tee temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        grep '[[:blank:]]-[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "-", $4, $5}' | ( [[ "$USEBP" ]] && tee -a temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals >> temp/"$fileID"fastaIntervals.bed
    elif [ "$INPUTBED" = "TRUE" ]; then
        grep '[[:blank:]]+[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "+", $7, $8}' | ( [[ "$USEBP" ]] && tee temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        grep '[[:blank:]]-[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "-", $7, $8}' | ( [[ "$USEBP" ]] && tee -a temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals >> temp/"$fileID"fastaIntervals.bed
    fi
    echo "Retrieving flanking FASTA sequence..."
    bedtools getfasta -fi $FASTAPATH -bed temp/"$fileID"fastaIntervals.bed -name -s > temp/"$fileID"seqToScan.FASTA
    if [ ! -s temp/"$fileID"seqToScan.FASTA ]; then
        echo "Error: no variants were returned following bedtools getfasta command. \n Exiting..."
        exit 1
    fi
    #seqScan: generates input strings for maxentscan and genesplicer as well as ESRseq scores
    echo "Scanning for motifs..."
    rm output/mesOmmitted/"$fileID" 2> /dev/null
    java -cp bin seqScan temp/"$fileID"seqToScan.FASTA -useESR $fileID 1>&2
    #run maxEntScan and confirm non-zero exit, since invalid inputs cause it to exit early
    echo "Running MaxEntScan..."
    perl score5.pl temp/"$fileID"mesDonorInput.txt | tee temp/"$fileID"mesDonorInputUnprocessed.txt | java -cp bin processScoresMES > temp/"$fileID"mesDonorScores.txt
    retVal=( ${PIPESTATUS[0]} )
    if [ $retVal -ne 0 ]; then
        echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
    exit $retVal
    fi
    perl score3.pl temp/"$fileID"mesAcceptorInput.txt | tee temp/"$fileID"mesAcceptorInputUnprocessed.txt | java -cp bin processScoresMES > temp/"$fileID"mesAcceptorScores.txt
    retVal=( ${PIPESTATUS[0]} )
    if [ $retVal -ne 0 ]; then
        echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
    exit $retVal
    fi
    #run genesplicer
    echo "Running GeneSplicer..."
    bin/linux/genesplicerAdapted temp/"$fileID"gsInput.FASTA human > temp/"$fileID"gsScores.txt
    #run branchpointer SNPs
    if [ "$USEBP" = "TRUE" -a "$USEBPINDELS" = "FALSE" ]; then
        echo "Running Branchpointer..."
        Rscript --slave --vanilla bin/bpProcessing.R "$fileID" $(pwd) "$ANNOTATION" "$GENOMEBUILD" &> output/"$fileID"bpLog.txt
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
        done < temp/"$fileID"bpInput.txt
        Rscript --slave --vanilla bin/bpProcessingINDELS.R "$fileID" $(pwd) "$ANNOTATION" "$GENOMEBUILD" &> output/"$fileID"bpIndelLog.txt
        #awk -v OFS=\\t '{print $2, $3, $4, $8, $9, $16, $17, $23, $24, $25, $26}' output/"$fileID"bpOutputIndels.txt > output/bpIndelsSummarised.txt
        #awk -v OFS=\\t '{print $2, $3, $4, $8, $9, $15, $16, $22, $23, $24, $25}' output/"$fileID"bpOutputSNPs.txt" > output/bpSNPsSummarised.txt
    fi
    #merge scores into one line
    echo "Processing scores..."
    cat temp/"$fileID"mesDonorScores.txt temp/"$fileID"mesAcceptorScores.txt temp/"$fileID"gsScores.txt temp/"$fileID"ESRoutput.txt data/"$gtfBasename"_SpliceSiteIntervals.txt sources/terminatingMergeLine.txt |
    sort -k1,1 -V -k 2,2n -k 3 -k 4 -s | tee teeTest2.txt | java -cp bin mergeOutput "$fileID"
    #sort predictions
    if [ -s temp/"$fileID"_donorCreating_unsorted.txt ]; then
        sort -gr -k11,11 temp/"$fileID"_donorCreating_unsorted.txt >> output/"$fileID"_donorCreating.txt
    else
        rm output/"$fileID"_donorCreating.txt
    fi 
    if [ -s temp/"$fileID"_acceptorCreating_unsorted.txt ]; then
        sort -gr -k11,11 temp/"$fileID"_acceptorCreating_unsorted.txt >> output/"$fileID"_acceptorCreating.txt
    else
        rm output/"$fileID"_acceptorCreating.txt
    fi 
    if [ -s temp/"$fileID"_withinSS_unsorted.txt ]; then
        sort -gr -k17,17 temp/"$fileID"_withinSS_unsorted.txt | cut -f1-15 >> output/"$fileID"_withinSS.txt
    else
        rm output/"$fileID"_withinSS.txt
    fi 
    #clean up temp files
    rm temp/"$fileID"* 2> /dev/null
done
