# Assess the effects of SNPs on branchpoint annotations (adapted from Branchpointer vignette)
# read in fileName as Rscript argument

args=(commandArgs(TRUE))
if(length(args)<4){
  print("Error: Not all arguments supplied.")
}else{
  querySNPFile <- paste(args[2], "/temp/", args[1], "bpInput.txt", sep='')
  outputPath <- paste(args[2], "/output/", args[1], sep='')
  annotationFile <- args[3]
  genomeBuild <- args[4]
  #fastaFile <- args[4]
  #bedtoolsLocation <- args[5]
  }
setwd(args[2])

#load libraries and exon annotation
suppressMessages(library(knitr))
suppressMessages(library(branchpointer))
if (identical(genomeBuild, "hg38")) {
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
} else
if (identical(genomeBuild, "hg19")) {
    suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
} else {
    print("Error: must specify genome build as hg38 or hg19")
}
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
#g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
exons <- gtfToExons(annotationFile)
#format and filter for branchpoint window SNPs
querySNP <- readQueryFile(querySNPFile, 
                          queryType = "SNP", 
                          exons = exons, 
                          filter = TRUE,
                          maxDist = 50)

#predict branchpoints and specify number of cores
branchpointPredictionsSNP <- predictBranchpoints(querySNP,
                                                 queryType = "SNP",
                                                 BSgenome = g,
                                                 useParallel = TRUE,
                                                 cores = 8)

#summarise and define prediction thresholds
querySNPSummary <- predictionsToSummary(querySNP,branchpointPredictionsSNP, probabilityCutoff = 0.52, probabilityChange = 0.15)

#identify branchpoints created/removed by SNP
if (FALSE) {
l <- length(querySNPSummary)
BPcreated_or_removed <- vector()
for (i in 1:l) {
  if (querySNPSummary$created_n[i]==1 | querySNPSummary$deleted_n[i]==1) {
    BPcreated_or_removed <- c(BPcreated_or_removed, i)
  }
}

#write to file
fileDestination <- paste(outputPath, "bpOutputSNPs.txt", sep='')
write.table(querySNPSummary[BPcreated_or_removed], fileDestination , quote = FALSE, sep="\t", row.names=FALSE)
}

fileDestination <- paste(outputPath, "bpOutputSNPs.txt", sep='')
write.table(querySNPSummary, fileDestination , quote = FALSE, sep="\t", row.names=FALSE)
