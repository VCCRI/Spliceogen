#!/bin/bash
function printVersion {
    printf "Spliceogen 2.0 1-October-2019\n"
}
function printHelp {
    cat <<-END
Usage:
------
3 required args:
1)  -input path/to/VCF/input/file(s).VCF
        Note: multiple input files are accepted "eg. -input *.vcf"
        Note: deprecated v1.0 tags "inputVCF" and "inputBED" are still accepted
2)  -gtf path/to/annotation.gtf
3)  -fasta path/to/genome.fa
optional arg:
4)  -branchpointer hgXX
        OR
    -branchpointerIndels hgXX
        Note: user must specify hg19 or hg38
END
}
#set default parameters
POSITIONAL=()
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
    -input)
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -inputVCF)
    INPUTFILES="$2"
    shift
    shift
    while [ "$1" ] && [[ ! $1 == *-* ]]; do
        INPUTFILES="$INPUTFILES $1"
        shift
    done
    ;;
    -inputBED)
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
#correct mismatches in "chr" nomenclature between gtf and fasta
gtfChr=$(zcat -f "$ANNOTATION" | grep -v '^GL000' | tail -1 | awk '{print $1}' | grep chr)
fastaChr=$(cat "$FASTAPATH" | head -1 | awk '{print $1}' | grep chr)
    gtfChrAdd=""
    gtfChrRemove="UnmatchedString"
    if [ "$fastaChr" != "" ]; then
        if [ "$gtfChr" == "" ]; then
            gtfChrAdd="chr"
        fi        
    elif [ "$fastaChr" == "" ]; then
        if [ "$gtfChr" != "" ]; then
            gtfChrRemove="chr"
        fi        
    fi
#prepare splice site intervals from annotation.gtf
gtfBasename=$(basename $ANNOTATION)
if [ ! -f data/"$gtfBasename"_SpliceSiteIntervals.txt ] || [[ "$ANNOTATION" -nt data/"$gtfBasename"_SpliceSiteIntervals.txt ]] ; then
    echo "Preparing splice site annotation..."
    zcat -f "$ANNOTATION" | grep '[[:blank:]]gene[[:blank:]]\|[[:blank:]]transcript[[:blank:]]\|[[:blank:]]exon[[:blank:]]' | grep -v '^GL000' | 
    sed "s/$gtfChrRemove//" | awk -v var="$gtfChrAdd" '{print var$0}' | java -cp bin getSpliceSiteIntervalsFromGTF > data/"$gtfBasename"_SpliceSiteIntervals.txt
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
    #check input file type
    FILETYPE=""
    nFields=$(zcat -f $FILE | tail -1 | wc -w)
    vcfHeader=$(zcat -f $FILE | head -1 | grep VCF)
    if [ "$nFields" -eq 4 ]; then
        FILETYPE="TSV"
    elif [ ! -z "$vcfHeader" ]; then
        FILETYPE="VCF"
    else
        FILETYPE="BED"
    fi
    #correct mismatches in "chr" nomenclature between variant input and provided gtf/fasta
    inputChr=$(zcat -f "$FILE" | tail -1 | awk '{print $1}' | grep chr)
    inputChrAdd=""
    inputChrRemove="UnmatchedString"
    if [ "$fastaChr" != "" ]; then
        if [ "$inputChr" == "" ]; then
            inputChrAdd="chr"
        fi        
    elif [ "$fastaChr" == "" ]; then
        if [ "$inputChr" != "" ]; then
            inputChrRemove="chr"
        fi        
    fi
    #sort body of input file
        zcat -f "$FILE" | grep "^#" > temp/"$fileID"_sorted
        if [ "$FILETYPE" == "TSV" ]; then
            zcat -f "$FILE" | grep -v "^#" | sort -k1,1 -k2,2n | sed "s/$inputChrRemove//" | awk -v OFS="\\t" -v var=$inputChrAdd '{print var$1, $2, $2, "x", "1", ".", $3, $4}' >> temp/"$fileID"_sorted
        else 
            zcat -f "$FILE" | grep -v "^#" | sort -k1,1 -k2,2n | sed "s/$inputChrRemove//" | awk -v OFS="\\t" -v var=$inputChrAdd '{print var$0}' >> temp/"$fileID"_sorted
        fi
    #check bedtools is installed
    bedtoolsLocation=$(which bedtools);
    if [ "$bedtoolsLocation" == "" ]; then
        printf -- 'Warning: Bedtools does not appear to be installed.\n';
        printf -- 'Get it here: https://bedtools.readthedocs.io/en/latest/content/installation.html\n';
    fi;
    #note: the bedtools getfasta "-name" behaviour is not backwards compatible. Between v2.26.0 and v2.27.0 the previous function  of "name" was given to "name+". So need to alter command for different versions
    recentBedtools=false
    bedtoolsVersion=$(bedtools -version)
    versionSort=$(echo -e "bedtools v2.26.0\nbedtools v2.27.0\n$bedtoolsVersion" | sort -V | head -2 | tail -1)
    if [ "$versionSort" == "bedtools v2.27.0" ]; then
        recentBedtools=true
    fi
    #bedtools intersect to get strand info
    echo "Retrieving strand info..."
    zcat -f "$ANNOTATION" | grep '[[:blank:]]gene[[:blank:]]' | sort -k1,1 -k4,4n | grep -v '^GL000' | awk -v var="$gtfChrAdd" -v OFS="\\t" '{print var$0}' | sed "s/$gtfChrRemove//" | bedtools intersect -a temp/"$fileID"_sorted -b stdin -wa -wb -sorted  > temp/"$fileID"unstrandedInput.txt 
    if [ $? -ne 0 ]; then
        echo "Warning. Bedtools intersect returned non-zero exit status. Intersection failed between provided variant VCF/BED file and provided GTF. See above error message for more details"
    fi
    if [ ! -s temp/"$fileID"unstrandedInput.txt ]; then
        echo "Error: no variants were returned following bedtools intersect between input file \""$fileID"\" and gtf. \n Exiting..."
        exit 1
    fi
    #generate flanking intervals.bed for bedtools getfasta and branchpointer input
    if [ "$FILETYPE" = "VCF" ]; then
        grep '[[:blank:]]+[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "+", $4, $5}' | ( [[ "$USEBP" ]] && tee temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        grep '[[:blank:]]-[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "-", $4, $5}' | ( [[ "$USEBP" ]] && tee -a temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals >> temp/"$fileID"fastaIntervals.bed
    elif [ "$FILETYPE" = "TSV" ] || [ "$FILETYPE" = "BED" ] ; then
        grep '[[:blank:]]+[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "+", $7, $8}' | ( [[ "$USEBP" ]] && tee temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals > temp/"$fileID"fastaIntervals.bed
        grep '[[:blank:]]-[[:blank:]]' temp/"$fileID"unstrandedInput.txt | awk -v OFS="\\t" '{print ".", $1, $2, "-", $7, $8}' | ( [[ "$USEBP" ]] && tee -a temp/"$fileID"bpInput.txt || cat ) | java -cp bin getFastaIntervals >> temp/"$fileID"fastaIntervals.bed
    fi
    echo "Retrieving flanking FASTA sequence..."
    if [ "$recentBedtools" == true ]; then
        bedtools getfasta -fi $FASTAPATH -bed temp/"$fileID"fastaIntervals.bed -name+ -s > temp/"$fileID"seqToScan.FASTA
    else
        bedtools getfasta -fi $FASTAPATH -bed temp/"$fileID"fastaIntervals.bed -name -s > temp/"$fileID"seqToScan.FASTA
    fi
    if [ ! -s temp/"$fileID"seqToScan.FASTA ]; then
        echo "Error: no variants were returned following bedtools getfasta command. \n Exiting..."
        exit 1
    fi
    #seqScan: generates input strings for maxentscan and genesplicer as well as ESRseq scores
    echo "Scanning for motifs..."
    rm output/"$fileID"mesOmmitted.txt 2> /dev/null
    rm output/"$fileID"refMismatch.txt 2> /dev/null
    java -cp bin seqScan temp/"$fileID"seqToScan.FASTA -useESR $fileID 1>&2
    if [ -s output/"$fileID"refMismatch.txt ]; then
        refMismatchCount=$(wc -l output/"$fileID"refMismatch.txt | awk '{print $1}')
        echo "Note: $refMismatchCount variants were excluded because the provided Reference allele does not match the nucleotide(s) in the provided FASTA. IDs of excluded variant(s) are outputted here: Spliceogen/output/""$fileID""refMismatch.txt"
    fi
    if [ -s output/"$fileID"mesOmmitted.txt ]; then
        mesOmmittedCount=$(wc -l output/"$fileID"mesOmmitted.txt | awk '{print $1}')
        echo "Note: $mesOmmittedCount variants were excluded from MaxEntScan because their flanking FASTA sequence contains invalid characters (most commonly \"n\"), which cannot be processed by MaxEntScan. IDs of ommitted variant(s) are listed in: Spliceogen/output/""$fileID""mesOmmitted.txt"
    fi
    #run maxEntScan and confirm non-zero exit, since invalid inputs cause it to exit early
    if [ -s temp/"$fileID"mesDonorInput.txt ] || [ -s temp/"$fileID"mesAcceptorInput.txt ] ; then
        echo "Running MaxEntScan..."
        perl score5.pl temp/"$fileID"mesDonorInput.txt | sed 's/;-;(-)/;/' | sed 's/;+;/;/' | sort -t '>' -k2 | java -cp bin processScoresMES > temp/"$fileID"mesDonorScores.txt
        retVal=( ${PIPESTATUS[0]} )
        if [ $retVal -ne 0 ]; then
            echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
        exit $retVal
        fi
        perl score3.pl temp/"$fileID"mesAcceptorInput.txt | sed 's/;-;(-)/;/' | sed 's/;+;/;/' | sort -t '>' -k2 | java -cp bin processScoresMES > temp/"$fileID"mesAcceptorScores.txt
        retVal=( ${PIPESTATUS[0]} )
        if [ $retVal -ne 0 ]; then
            echo "MaxEntScan returned non-zero exit status. It is likely not all variants were processed. Exiting..."
        exit $retVal
        fi
    else
        echo "No input for MaxEntScan"
    fi
    #run genesplicer
    if [ -s temp/"$fileID"gsInput.FASTA ] ; then
        echo "Running GeneSplicer..."
        bin/linux/genesplicerAdapted temp/"$fileID"gsInput.FASTA human > temp/"$fileID"gsScores.txt
    else
        echo "No input for GeneSplicer"
    fi
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
    scoresToMerge=""
    if [ -s temp/"$fileID"gsScores.txt ] ; then
        scoresToMerge="temp/"$fileID"gsScores.txt"
    fi
    if [ -s temp/"$fileID"mesDonorScores.txt ] ; then
        scoresToMerge="$scoresToMerge temp/"$fileID"mesDonorScores.txt"
    fi
    if [ -s temp/"$fileID"mesAcceptorScores.txt ] ; then
        scoresToMerge="$scoresToMerge temp/"$fileID"mesAcceptorScores.txt"
    fi
    if [ -s temp/"$fileID"ESRoutput.txt ] ; then
        scoresToMerge="$scoresToMerge temp/"$fileID"ESRoutput.txt"
    fi
    #edit splice site intervals file in event of changed input/fasta "chr" nomenclature
    intervalsFileChr=$(cat data/"$gtfBasename"_SpliceSiteIntervals.txt | head -1 | awk '{print $1}' | grep chr)
    if [ "$fastaChr" != "" ]; then
        if [ "$intervalsFileChr" == "" ]; then
            sed -i 's/^/chr/' data/"$gtfBasename"_SpliceSiteIntervals.txt
        fi
    elif [ "$fastaChr" == "" ]; then
        if [ "$gtfChr" != "" ]; then
            sed -i "s/^chr//" data/"$gtfBasename"_SpliceSiteIntervals.txt
        fi
    fi
    checkScoresExist=$(echo "$scoresToMerge" | grep "temp")
    if [ -z "$checkScoresExist" ]; then
        echo "No MaxEntScan/GeneSplicer/ESRseq scores to process"
    else
        echo "Processing scores..."
        cat $(echo "$scoresToMerge") data/"$gtfBasename"_SpliceSiteIntervals.txt sources/terminatingMergeLine.txt | sort -k1,1 -V -k 2,2n -k 3 -k 4 -s | java -cp bin mergeOutput "$fileID" inputAdd="$inputChrRemove" inputRemove="$inputChrAdd"
    fi
    #sort predictions
    if [ -s temp/"$fileID"_gain_unsorted.txt ]; then
        echo -e "#CHR\tSTART\tREF\tALT\tGENE\tdonGainP\taccGainP" > output/"$fileID"_ssGain.txt
        sort -gr -k8,8 temp/"$fileID"_gain_unsorted.txt | cut -f1-7 >> output/"$fileID"_ssGain.txt
    fi 
    if [ -s temp/"$fileID"_loss_unsorted.txt ]; then
        echo -e "#CHR\tSTART\tREF\tALT\tGENE\twithinSS\tdonLossP\taccLossP" > output/"$fileID"_withinSS.txt
        sort -gr -k9,9 temp/"$fileID"_loss_unsorted.txt | cut -f1-8 >> output/"$fileID"_withinSS.txt
    fi 
    #clean up temp files
    rm temp/"$fileID"* 2> /dev/null
done
