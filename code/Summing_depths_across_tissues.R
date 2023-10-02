#!/usr/bin/env R
#Summing_depths_across_tissues.R



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
# mat_outputFile: Path to which the summed counts across tissues for haplotype 1
#     should be written.
# pat_outputFile: Path to which the summed counts across tissues for haplotype 2
#     should be written.
#
#
################################################################################


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################


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
mat_outputFile <- args[3]
pat_outputFile <- args[4]

#matDepths <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg.txt"
#patDepths <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg.txt"
#mat_outputFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#pat_outputFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"



matDepths <- read.table(matDepths, header=T, sep="\t", stringsAsFactors=F)
patDepths <- read.table(patDepths, header=T, sep="\t", stringsAsFactors=F)

colnames(matDepths) <- toupper(colnames(matDepths))
colnames(patDepths) <- toupper(colnames(patDepths))


allExprs <- colnames(matDepths)
allExprs <- matrix(nrow=length(allExprs), ncol=3, data=unlist(strsplit(allExprs, split="_")), byrow=T)

allTFs <- unique(allExprs[,2])


summed_matTable <- matrix(nrow=nrow(matDepths), ncol=length(allTFs), data=1)
summed_patTable <- matrix(nrow=nrow(patDepths), ncol=length(allTFs), data=1)
colnames(summed_matTable) <- allTFs
rownames(summed_matTable) <- rownames(matDepths)
colnames(summed_patTable) <- allTFs
rownames(summed_patTable) <- rownames(patDepths)

for (i in 1:length(allTFs)) {
  print(paste(i, length(allTFs)))
  if(length(grep(paste(allTFs[i], "_", sep=""), colnames(matDepths)))>1) {
    matSums <- rowSums(matDepths[,grep(paste(allTFs[i], "_", sep=""), colnames(matDepths))])
    patSums <- rowSums(patDepths[,grep(paste(allTFs[i], "_", sep=""), colnames(patDepths))])
  }
  if(length(grep(paste(allTFs[i], "_", sep=""), colnames(matDepths)))==1) {
    matSums <- matDepths[,grep(paste(allTFs[i], "_", sep=""), colnames(matDepths))]
    patSums <- patDepths[,grep(paste(allTFs[i], "_", sep=""), colnames(patDepths))]
  }
  summed_matTable[,i] <- matSums
  summed_patTable[,i] <- patSums

}

write.table(summed_matTable, mat_outputFile, row.names=T, col.names=T, sep="\t", quote=F)
write.table(summed_patTable, pat_outputFile, row.names=T, col.names=T, sep="\t", quote=F)







#
