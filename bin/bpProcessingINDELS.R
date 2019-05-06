# Assess the effects of SNPs and INDELS on branchpoint annotations (adapted from Branchpointer vignette)
# read in fileName as Rscript argument
args=(commandArgs(TRUE))
if(length(args)<4){
  print("Error: Not all arguments supplied.")
}else{
  querySNPFile <- paste(args[2], "/temp/", args[1], "bpInput.txt", sep='')
  queryIndelFile <- paste(args[2], "/temp/", args[1], "bpInputIndels.txt", sep='')
  outputPath <- paste(args[2], "/output/", args[1], sep='')
  annotationFile <- args[3]
  genomeBuild <- args[4]
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
exons <- gtfToExons(annotationFile)

#handle indels (if file not empty)
if (file.exists(queryIndelFile)) {
  #format and filter for branchpoint window indels
  queryIndel <- readQueryFile(queryIndelFile, 
                          queryType = "indel", 
                          exons = exons)

  #predict branchpoints and specify number of cores
  branchpointPredictionsIndel <- predictBranchpoints(queryIndel,
                                                 queryType = "indel",
                                                 BSgenome = g,
                                                 useParallel = TRUE,
                                                 cores = 8)

  #summarise and define prediction thresholds
  queryIndelSummary <- predictionsToSummary(queryIndel,branchpointPredictionsIndel, probabilityCutoff = 0.52, probabilityChange = 0.15)
  #write to file
  fileDestinationIndels <- paste(outputPath, "bpOutputIndels.txt", sep='')
  write.table(queryIndelSummary, fileDestinationIndels , quote = FALSE, sep="\t", row.names=FALSE)
}
#handle SNPs
if (file.exists(querySNPFile)) {
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
  #write to file
  fileDestinationSNPs <- paste(outputPath, "bpOutputSNPs.txt", sep='')
  write.table(querySNPSummary, fileDestinationSNPs , quote = FALSE, sep="\t", row.names=FALSE)
}

#optional: filter only for branchpoints created/removed by Indel
if (FALSE) {
  l <- length(queryIndelSummary)
BPcreated_or_removed <- vector()
for (i in 1:l) {
  if (queryIndelSummary$created_n[i]>=1 | queryIndelSummary$deleted_n[i]>=1) {
    BPcreated_or_removed <- c(BPcreated_or_removed, i)
  }
}
}
