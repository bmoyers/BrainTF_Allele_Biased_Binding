#!/usr/bin/env R
#CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R


################################################################################
#This script is used to compile all p-value quantifications over all TFs when
#     different datasets for the same tissue are summed across tissues. Note that
#     it is expected that the column names are formatted as "TISSUE_TF" for
#     proper processing.
#
#Run under R 4.1.0, but can be used under other versions.
#
#This program takes the following arguments:
# matDepths: Path to the full counts table of region by experiment for haplotype 1,
#     produced by the script CompilingHaplotypeReadsTable.R
# patDepths: Path to the full counts table of region by experiment for haplotype 2,
#     produced by the script CompilingHaplotypeReadsTable.R
# outputFile: Path to which the final summed-across-tissues pvalue table should be
#     written.
#
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
#Define Functions
################################################################################
################################################################################



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################


args <- commandArgs(trailingOnly=T)
matDepths <- args[1]
patDepths <- args[2]
outputFile <- args[3]


matDepths <- read.table(matDepths, header=T, sep="\t", stringsAsFactors=F)
patDepths <- read.table(patDepths, header=T, sep="\t", stringsAsFactors=F)




colnames(matDepths) <- toupper(colnames(matDepths))
colnames(patDepths) <- toupper(colnames(patDepths))


allTFs <- colnames(matDepths)


summedPvalTable <- matrix(nrow=nrow(matDepths), ncol=length(allTFs), data=1)
colnames(summedPvalTable) <- allTFs
rownames(summedPvalTable) <- rownames(matDepths)
for (i in 1:length(allTFs)) {
  print(paste(i, length(allTFs)))
  matSums <- matDepths[,allTFs[i]]
  patSums <- patDepths[,allTFs[i]]


  toReplace <- which((matSums+patSums)>=6)
  theDepths <- cbind(rowSums(cbind(matSums, patSums)))
  test_vector <- floor(apply(cbind(matSums, patSums), 1, max))


  bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
  thePvals_replace <- mapply(bt, test_vector[toReplace], theDepths[toReplace,1])
  thePvals_replace <- as.numeric(thePvals_replace)
  summedPvalTable[toReplace,i] <- thePvals_replace

}


write.table(summedPvalTable, outputFile, row.names=T, col.names=T, sep="\t", quote=F)
