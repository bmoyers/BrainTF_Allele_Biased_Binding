#!/usr/bin/env R
#CompilingHaplotypePvalsTable.R

################################################################################
#This script is used to compile all p-value quantifications over all experiments
#     into a single table for ease of reference and use by later scripts.
#
#Run under R 4.1.0, but can be used under other versions.
#
#This program takes the following arguments:
# theDir: Path to the directory storying the p-value tables for individual experiments,
#     produced by the DonorSpecificMapping.sh script.
# pvalOutput: Path to which the final p-value table will be written.
#
#This program has the following optional arguments:
# toRemove: All arguments after "pvalOutput" should be in the form of
#     experiment names, e.g. "AMY_BCL11A", and will be removed from consideration.
#     This was done because certain experiments were flagged by collaborators
#     as problematic after I had performed mapping.
#
################################################################################

################################################################################
################################################################################
#Load Libraries and define options.
################################################################################
################################################################################

options("scipen"=100)


################################################################################
################################################################################
#Define Functions.
################################################################################
################################################################################



################################################################################
################################################################################
#Begin Script.
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
theDir <- args[1]
pvalOutput <- args[2]
if(length(args)>2) {
  toRemove <- args[3:length(args)]
}
if(length(args)==2) {
  toRemove <- NA
}



thePvalFiles <- list.files(theDir, pattern="_binomial_pval.txt", full.names=T)

firstPvals <- read.table(thePvalFiles[1])

finalPvals <- matrix(nrow=nrow(firstPvals), ncol=length(thePvalFiles), data=0)

finalPvals[,1] <- firstPvals[,1]

pval_colnames <- c(colnames(firstPvals)[1])

for (i in 2:length(thePvalFiles)) {
  print(paste(i, length(thePvalFiles), sep="  "))

  thisPval <- read.table(thePvalFiles[i])
  finalPvals[,i] <- thisPval[,1]
  pval_colnames <- c(pval_colnames, colnames(thisPval)[1])

}

colnames(finalPvals) <- pval_colnames

rownames(finalPvals) <- rownames(firstPvals)

colnames(finalPvals) <- toupper(colnames(finalPvals))


if(!is.na(toRemove)) {
  for (i in 1:length(toRemove)) {
    if(length(grep(toRemove[i], colnames(finalPvals)))>0) { finalPvals <- finalPvals[,-grep(toRemove[i], colnames(finalPvals))]}

  }
}


write.table(finalPvals, pvalOutput, sep="\t", quote=F, row.names=T, col.names=T)
