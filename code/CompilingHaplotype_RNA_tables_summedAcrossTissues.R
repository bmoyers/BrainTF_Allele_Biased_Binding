#CompilingHaplotype_RNA_tables_summedAcrossTissues.R


################################################################################
#This script is used to compile all p-value quantifications over all TFs when
#     different datasets for the same tissue are summed across tissues. Note that
#     it is expected that the column names are formatted as "TISSUE_TF" for
#     proper processing.
#
#Run under R 4.1.0, but can be used under other versions.
#
#This program takes the following arguments:
# outDir: Directory to which the final summed-across-tissues pvalue table  and
#     depth tables should be written.
# rna_dir: Directory to which the RNA depth tables were written.
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

outDir <- args[1]
rna_dir <- args[2]

#outDir <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/"
#rna_dir <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/RNA_Mapped_vg"




rnaFiles <- list.files(rna_dir, full.names=T)
matFiles <- rnaFiles[grep("maternal", rnaFiles)]
patFiles <- rnaFiles[grep("paternal", rnaFiles)]



allTissues <- c()
for (i in 1:length(matFiles)) {
  thisFile <- strsplit(matFiles[i], split="/")[[1]]
  thisFile <- thisFile[length(thisFile)]
  thisTissue <- strsplit(thisFile, split="_")[[1]][2]
  allTissues <- c(allTissues, thisTissue)
}
allTissues <- unique(allTissues)

testMat <- read.table(matFiles[1], header=T, sep="\t", stringsAsFactors=F)
finalMatCounts <- matrix(nrow=nrow(testMat), ncol=length(allTissues), data=0)
colnames(finalMatCounts) <- toupper(allTissues)
rownames(finalMatCounts) <- rownames(testMat)
finalPatCounts <- finalMatCounts
finalPvals <- matrix(nrow=nrow(testMat), ncol=length(allTissues), data=1)
colnames(finalPvals) <- toupper(allTissues)
rownames(finalPvals) <- rownames(testMat)

for (i in 1:length(allTissues)) {
  theseMats <- matFiles[grep(paste("_", allTissues[i], "_", sep=""), matFiles)]
  matCounts <- read.table(theseMats[1], header=T, sep="\t", stringsAsFactors=F)
  for (j in 2:length(theseMats)) {
    this_matCounts <- read.table(theseMats[j], header=T, sep="\t", stringsAsFactors=F)
    matCounts <- matCounts+this_matCounts
  }
  thesePats <- patFiles[grep(paste("_", allTissues[i], "_", sep=""), patFiles)]
  patCounts <- read.table(thesePats[1], header=T, sep="\t", stringsAsFactors=F)
  for (j in 2:length(thesePats)) {
    this_patCounts <- read.table(thesePats[j], header=T, sep="\t", stringsAsFactors=F)
    patCounts <- patCounts+this_patCounts
  }
  finalMatCounts[,i] <- matCounts[,1]
  finalPatCounts[,i] <- patCounts[,1]


  toReplace <- which((matCounts[,1]+patCounts[,1])>=6)
  theDepths <- cbind(rowSums(cbind(matCounts, patCounts)))
  test_vector <- floor(apply(cbind(matCounts[,1], patCounts[,1]), 1, max))

  bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
  thePvals_replace <- mapply(bt, test_vector[toReplace], theDepths[toReplace,1])
  thePvals_replace <- as.numeric(thePvals_replace)
  finalPvals[toReplace,i] <- thePvals_replace




}

pvals_summed <- rep(1, nrow(finalPvals))
matSums <- rowSums(finalMatCounts)
patSums <- rowSums(finalPatCounts)

toReplace <- which((matSums+patSums)>=6)
theDepths <- cbind(rowSums(cbind(matSums, patSums)))
test_vector <- floor(apply(cbind(matSums, patSums), 1, max))


bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
thePvals_replace <- mapply(bt, test_vector[toReplace], theDepths[toReplace,1])
thePvals_replace <- as.numeric(thePvals_replace)
pvals_summed[toReplace] <- thePvals_replace

finalPvals <- cbind(finalPvals, pvals_summed)
finalMatCounts <- cbind(finalMatCounts, matSums)
finalPatCounts <- cbind(finalPatCounts, patSums)

write.table(finalMatCounts, paste(outDir, "HaplotypeReadsTable_RNA_maternal_vg.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
write.table(finalPatCounts, paste(outDir, "HaplotypeReadsTable_RNA_paternal_vg.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
write.table(finalPvals, paste(outDir, "HaplotypePvalsTable_RNA_binom_vg.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F)
