# Assess the effects of SNPs and INDELS on branchpoint annotations (adapted from Branchpointer vignette)
# read in fileName as Rscript argument
args=(commandArgs(TRUE))
if(length(args)<3){
  print("Error: Not all arguments supplied.")
}else{
  querySNPFile <- args[1]
  queryIndelFile <- args[2]
  outputPath <- args[3]
  annotationPath <- args[4]
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

#identify branchpoints created/removed by Indel
if (FALSE) {
  l <- length(queryIndelSummary)
BPcreated_or_removed <- vector()
for (i in 1:l) {
  if (queryIndelSummary$created_n[i]>=1 | queryIndelSummary$deleted_n[i]>=1) {
    BPcreated_or_removed <- c(BPcreated_or_removed, i)
  }
}
}

#write to file
setwd(outputPath)
outputFileSNPs <- "./output/bpOutputSNPs.txt"
outputFileIndels <- "./output/bpOutputIndels.txt"
write.table(querySNPSummary, outputFileSNPs , quote = FALSE, sep="\t")
write.table(queryIndelSummary, outputFileIndels , quote = FALSE, sep="\t")
