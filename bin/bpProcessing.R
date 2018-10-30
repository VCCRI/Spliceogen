# Assess the effects of SNPs on branchpoint annotations (adapted from Branchpointer vignette)
# read in fileName as Rscript argument
args=(commandArgs(TRUE))
if(length(args)<3){
  print("Error: Not all arguments supplied.")
}else{
  querySNPFile <- args[1]
  outputPath <- args[2]
  annotationPath <- args[3]
  }

#load libraries and exon annotation
suppressMessages(library(knitr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(branchpointer))
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
exons <- gtfToExons(annotationPath)
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
fileDestination <- "bpOutputSNPs.txt"
write.table(querySNPSummary[BPcreated_or_removed], fileDestination , quote = FALSE, sep="\t")
}

setwd(outputPath)
fileDestination <- "./output/bpOutputSNPs.txt"
write.table(querySNPSummary, fileDestination , quote = FALSE, sep="\t")
