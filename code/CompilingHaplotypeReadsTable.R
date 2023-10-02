#!/usr/bin/env R
#CompilingHaplotypeReadsTable.R


################################################################################
#This script is used to compile all count quantifications over all experiments
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





#args <- commandArgs(trailingOnly=T)
#theDir <- args[1]
#matOutput <- args[2]
#patOutput <- args[3]



################################################################################
#Donor1
################################################################################



theDir <- "/gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/MappedReads_vg"
matOutput <- "/gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg.txt"
patOutput <- "/gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg.txt"

thePatFiles <- list.files(theDir, pattern="counts_paternal.txt", full.names=T)
theMatFiles <- list.files(theDir, pattern="counts_maternal.txt", full.names=T)

firstPatCounts <- read.table(thePatFiles[1])
firstMatCounts <- read.table(theMatFiles[1])

finalPatCounts <- matrix(nrow=nrow(firstPatCounts), ncol=length(thePatFiles), data=0)
finalMatCounts <- matrix(nrow=nrow(firstMatCounts), ncol=length(theMatFiles), data=0)

finalPatCounts[,1] <- firstPatCounts[,1]
finalMatCounts[,1] <- firstMatCounts[,1]

pat_colnames <- c(gsub("_variantReads_counts", "", colnames(firstPatCounts)[1]))
mat_colnames <- c(gsub("_variantReads_counts", "", colnames(firstMatCounts)[1]))

for (i in 2:length(thePatFiles)) {
  print(paste(i, length(thePatFiles), sep="  "))

  thisPat <- read.table(thePatFiles[i])
  finalPatCounts[,i] <- thisPat[,1]
  thisPatColnames <- gsub("_variantReads_counts", "", colnames(thisPat)[1])
  if(thisPatColnames%in%colnames(pat_colnames)) {
    print(i)
    break()
  }
  pat_colnames <- c(pat_colnames, thisPatColnames)

  thisMat <- read.table(theMatFiles[i])
  finalMatCounts[,i] <- thisMat[,1]
  mat_colnames <- c(mat_colnames, gsub("_variantReads_counts", "", colnames(thisMat)[1]))

}

finalPatCounts_safe <- finalPatCounts
finalMatCounts_safe <- finalMatCounts

colnames(finalPatCounts) <- pat_colnames
colnames(finalMatCounts) <- mat_colnames

rownames(finalPatCounts) <- rownames(firstPatCounts)
rownames(finalMatCounts) <- rownames(firstMatCounts)



allMat <- matrix(nrow=ncol(finalPatCounts), ncol=3, data=unlist(strsplit(colnames(finalPatCounts), split="_")), byrow=T)
allTFs <- unique(allMat[,2])
allTissues <- unique(allMat[,1])
a <- table(allMat[,2])
a <- a[order(as.numeric(a), decreasing=T)]



colnames(finalPatCounts) <- toupper(colnames(finalPatCounts))
colnames(finalMatCounts) <- toupper(colnames(finalMatCounts))


#OLIG2 1224 (donor 1) were ALL marked to not use
if(length(grep("OLIG2", colnames(finalPatCounts)))>0) { finalPatCounts <- finalPatCounts[,-grep("OLIG2", colnames(finalPatCounts))]}
if(length(grep("OLIG2", colnames(finalMatCounts)))>0) { finalMatCounts <- finalMatCounts[,-grep("OLIG2", colnames(finalMatCounts))]}


#ASCL1_CB 1224 (donor 1) is marked to not use
#colnames(finalPatCounts)[grep("CB_ASCL1", colnames(finalPatCounts))]
if(length(grep("CB_ASCL1", colnames(finalPatCounts)))>0) { finalPatCounts <- finalPatCounts[,-grep("CB_ASCL1", colnames(finalPatCounts))]}
if(length(grep("CB_ASCL1", colnames(finalMatCounts)))>0) { finalMatCounts <- finalMatCounts[,-grep("CB_ASCL1", colnames(finalMatCounts))]}

#ATF7_CB 1224 (donor1) is marked to not use
colnames(finalPatCounts)[grep("CB_ATF7", colnames(finalPatCounts))]
if(length(grep("CB_ATF7", colnames(finalPatCounts)))>0) { finalPatCounts <- finalPatCounts[,-grep("CB_ATF7", colnames(finalPatCounts))]}
if(length(grep("CB_ATF7", colnames(finalMatCounts)))>0) { finalMatCounts <- finalMatCounts[,-grep("CB_ATF7", colnames(finalMatCounts))]}




write.table(finalPatCounts, patOutput, sep="\t", quote=F, row.names=T, col.names=T)
write.table(finalMatCounts, matOutput, sep="\t", quote=F, row.names=T, col.names=T)






################################################################################
#Donor2
################################################################################






theDir <- "/gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/MappedReads_vg"
matOutput <- "/gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg.txt"
patOutput <- "/gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg.txt"

thePatFiles <- list.files(theDir, pattern="counts_paternal.txt", full.names=T)
theMatFiles <- list.files(theDir, pattern="counts_maternal.txt", full.names=T)

firstPatCounts <- read.table(thePatFiles[1])
firstMatCounts <- read.table(theMatFiles[1])

finalPatCounts <- matrix(nrow=nrow(firstPatCounts), ncol=length(thePatFiles), data=0)
finalMatCounts <- matrix(nrow=nrow(firstMatCounts), ncol=length(theMatFiles), data=0)

finalPatCounts[,1] <- firstPatCounts[,1]
finalMatCounts[,1] <- firstMatCounts[,1]

pat_colnames <- c(gsub("_variantReads_counts", "", colnames(firstPatCounts)[1]))
mat_colnames <- c(gsub("_variantReads_counts", "", colnames(firstMatCounts)[1]))

for (i in 2:length(thePatFiles)) {
  print(paste(i, length(thePatFiles), sep="  "))

  thisPat <- read.table(thePatFiles[i])
  finalPatCounts[,i] <- thisPat[,1]
  pat_colnames <- c(pat_colnames, gsub("_variantReads_counts", "", colnames(thisPat)[1]))

  thisMat <- read.table(theMatFiles[i])
  finalMatCounts[,i] <- thisMat[,1]
  mat_colnames <- c(mat_colnames, gsub("_variantReads_counts", "", colnames(thisMat)[1]))

}


finalPatCounts_safe <- finalPatCounts
finalMatCounts_safe <- finalMatCounts

colnames(finalPatCounts) <- pat_colnames
colnames(finalMatCounts) <- mat_colnames

rownames(finalPatCounts) <- rownames(firstPatCounts)
rownames(finalMatCounts) <- rownames(firstMatCounts)




allMat <- matrix(nrow=ncol(finalPatCounts), ncol=3, data=unlist(strsplit(colnames(finalPatCounts), split="_")), byrow=T)
allTFs <- unique(allMat[,2])
allTissues <- unique(allMat[,1])
a <- table(allMat[,2])
a <- a[order(as.numeric(a), decreasing=T)]


colnames(finalPatCounts) <- toupper(colnames(finalPatCounts))
colnames(finalMatCounts) <- toupper(colnames(finalMatCounts))





write.table(finalPatCounts, patOutput, sep="\t", quote=F, row.names=T, col.names=T)
write.table(finalMatCounts, matOutput, sep="\t", quote=F, row.names=T, col.names=T)
