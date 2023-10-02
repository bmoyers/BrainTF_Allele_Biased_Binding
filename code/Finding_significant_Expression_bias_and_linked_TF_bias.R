#!/usr/bin/env R
#Finding_significant_Expression_bias_and_linked_TF_bias.R



################################################################################
#This script is used to identify cases of RNA bias which are in the same phase
#     as a variant with TF association bias. It then saves a table of these
#     variant pairs.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges, and matrixStats packages are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# rna_sig_d1: A table containing information about RNA-biased alleles in donor 1, produced
#     by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R
# rna_mat_d1: A table containing the depth of coverge for Haplotype 1 in donor 1
#     for RNA reads, produced by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R
# rna_pat_d1: A table containing the depth of coverge for Haplotype 2 in donor 1
#     for RNA reads, produced by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R
# phase_grouping_d1: the vcf for donor 1. Provided via the following DOI:
#     https://doi.org/10.7303/syn4921369
# tf_sig_d1: A table containing information about TF-biased alleles summed across tissues for donor 1,
#     produced by the script CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# tf_mat_d1: A table containing the depth of coverge for Haplotype 1 in donor 1
#     for TF reads summed across tissues, produced by the script
#     Summing_depths_across_tissues.R
# tf_pat_d1: A table containing the depth of coverge for Haplotype 2 in donor 1
#     for TF reads summed across tissues, produced by the script
#     Summing_depths_across_tissues.R
# rna_sig_d2: A table containing information about RNA-biased alleles in donor 2, produced
#     by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R
# rna_mat_d2: A table containing the depth of coverge for Haplotype 1 in donor 2
#     for RNA reads, produced by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R
# rna_pat_d2: A table containing the depth of coverge for Haplotype 2 in donor 2
#     for RNA reads, produced by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R
# phase_grouping_d2: the vcf for donor 2, Provided via the following DOI:
#     https://doi.org/10.7303/syn4921369
# tf_sig_d2: A table containing information about TF-biased alleles summed across tissues for donor 2,
#     produced by the script CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# tf_mat_d2: A table containing the depth of coverge for Haplotype 1 in donor 2
#     for TF reads summed across tissues, produced by the script
#     Summing_depths_across_tissues.R
# tf_pat_d2: A table containing the depth of coverge for Haplotype 2 in donor 2
#     for TF reads summed across tissues, produced by the script
#     Summing_depths_across_tissues.R
#
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################



library(GenomicRanges)
library(matrixStats)


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

#This function determines the phase group for each of the variants
#within the table.
getPhaseGroup <- function(sig_rna_regions, phase_grouping) {

  phase_grouping[,3] <- phase_grouping[,2]
  colnames(phase_grouping) <- c("chr", "start", "end", "a1", "a2", "phase_group", "haplotype")
  phase_grouping <- phase_grouping[phase_grouping[,"phase_group"]!=".",]
  phase_grouping <- as.data.frame(phase_grouping)
  phase_grouping[,1] <- as.character(phase_grouping[,1])
  phase_grouping[,2] <- as.numeric(as.character(phase_grouping[,2]))
  phase_grouping[,3] <- as.numeric(as.character(phase_grouping[,3]))
  phase_grouping[,4] <- as.character(phase_grouping[,4])
  phase_grouping[,5] <- as.character(phase_grouping[,5])
  phase_grouping[,6] <- as.numeric(as.character(phase_grouping[,6]))
  phase_grouping[,7] <- as.character(phase_grouping[,7])

  gr_phase <- makeGRangesFromDataFrame(phase_grouping, keep.extra.columns=TRUE)

  sig_rna_regions <- as.data.frame(sig_rna_regions)
  sig_rna_regions[,1] <- as.character(sig_rna_regions[,1])
  sig_rna_regions[,2] <- as.numeric(as.character(sig_rna_regions[,2]))
  sig_rna_regions[,3] <- as.numeric(as.character(sig_rna_regions[,3]))
  colnames(sig_rna_regions) <- c("chr", "start", "end")

  gr_sig_rna <- makeGRangesFromDataFrame(sig_rna_regions, keep.extra.columns=TRUE)



  theVars <- findOverlaps(gr_phase, gr_sig_rna, type="within", select="all", ignore.strand=T)
  theVars <- as.data.frame(theVars)
  theVars <- theVars[order(theVars[,2]),]
  finalPhases <- phase_grouping[theVars[,1],"phase_group"]
  theVars <- cbind(theVars, finalPhases)

  phases <- c()
  for (i in 1:nrow(sig_rna_regions)) {
    thesePhases <- unique(theVars[theVars[,2]==i,3])
    if(length(thesePhases)==1) {phases <- c(phases, thesePhases)}
    if(length(thesePhases)>1) {
      print(i)
      break()
    }
  }

  sig_rna_regions <- cbind(sig_rna_regions, phases)

  return(sig_rna_regions)

}




#This function combines cases of RNA and TF bias which are in different locations
#but within the same phase.
getPhasedTFCases <- function(rna_sig_regions_wPhase, tf_sig_regions_wPhase, saveFile) {
  #rna_sig_regions_wPhase <- rna_sig_d2_regions_wPhase
  #tf_sig_regions_wPhase <- tf_sig_d2_regions_wPhase
  #saveFile <- paste(outDir, "rna_sig_d2_regions_wPhase_w_SignificantTFInPhase_alt.txt", sep="")

  theChroms <- unique(rna_sig_regions_wPhase[,1])

  #retMat <- c()
  for (i in 1:length(theChroms)) {
    print(paste(theChroms[i], i, length(theChroms)))
    this_rna_sig_regions_wPhase <- rna_sig_regions_wPhase[rna_sig_regions_wPhase[,1]==theChroms[i],]
    this_tf_sig_regions_wPhase <- tf_sig_regions_wPhase[tf_sig_regions_wPhase[,1]==theChroms[i],]
    uniquePhases <- unique(this_rna_sig_regions_wPhase[,4])
    uniquePhases <- uniquePhases[uniquePhases%in%this_tf_sig_regions_wPhase[,4]]
    this_rna_sig_regions_wPhase <- this_rna_sig_regions_wPhase[this_rna_sig_regions_wPhase$phase%in%uniquePhases,]
    this_tf_sig_regions_wPhase <- this_tf_sig_regions_wPhase[this_tf_sig_regions_wPhase$phase%in%uniquePhases,]

    a <- lapply(this_rna_sig_regions_wPhase$phases, function(x) which (this_tf_sig_regions_wPhase$phases==x))
    matchMat <- matrix(nrow=length(unlist(a)), ncol=2, data=NA)
    matchMat[,2] <- unlist(a)
    currentRow <-1
    for (j in 1:length(a)) {
      thisSet <- a[[j]]
      thisSet <- rep(j, length(thisSet))
      endRow <- currentRow+length(thisSet)-1
      matchMat[currentRow:endRow,1] <- thisSet
      currentRow <- endRow+1
    }
    thisSet <- cbind(this_rna_sig_regions_wPhase[matchMat[,1],], this_tf_sig_regions_wPhase[matchMat[,2],])
    if(nrow(thisSet[thisSet[,4]==thisSet[,8],])!=nrow(thisSet)) {
      print(paste("Error with phase matching ", theChroms[i]))
      break()
    }
    thisSet_lines <- paste(thisSet[,1], thisSet[,2], thisSet[,3], thisSet[,5], thisSet[,6], thisSet[,7], thisSet[,8], sep="\t")
    write(thisSet_lines, saveFile, append=T, sep="\n")
  }
  retMat <- read.table(saveFile, header=F, sep="\t", stringsAsFactors=F)
  colnames(retMat) <- c("chr", "start", "end", "chr_tf", "start_tf", "end_tf", "phase")
  system(paste("rm ", saveFile, sep=""))

  return(retMat)

}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################


outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/"

rna_sig_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_RNA_binom_vg.txt"
rna_mat_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_RNA_maternal_vg.txt"
rna_pat_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_RNA_paternal_vg.txt"
phase_grouping_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants.vcf.table"
tf_sig_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
tf_mat_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
tf_pat_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"


rna_sig_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_RNA_binom_vg.txt"
rna_mat_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_RNA_maternal_vg.txt"
rna_pat_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_RNA_paternal_vg.txt"
phase_grouping_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.table"
tf_sig_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
tf_mat_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
tf_pat_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"



thePromoters <- read.table(thePromoters, header=F, sep="\t", stringsAsFactors=F)


###########
#For each significant RNA region:
#1) Determine whether it's Maternal or Paternally preferred.
#2) Get the Phase Grouping.
#3) Identify cases where there is a Phased variant that is significant
#   for a TF.
#3b) Determine the TF preference for Mat or Pat
#4) Identify the closest TSS (and distance)
$5) Save.
###########

################################################################################
################################################################################
#Donor 1
################################################################################
################################################################################

################################################################################
#Read in the RNA tables. Remove X chrom. Restrict to significant.
################################################################################

rna_sig_d1 <- read.table(rna_sig_d1, header=T, sep="\t", stringsAsFactors=F)
rna_mat_d1 <- read.table(rna_mat_d1, header=T, sep="\t", stringsAsFactors=F)
rna_pat_d1 <- read.table(rna_pat_d1, header=T, sep="\t", stringsAsFactors=F)


toRemove <- grep("chrX", rownames(rna_sig_d1))
if(length(toRemove)>0) {
  rna_sig_d1 <- rna_sig_d1[-grep("chrX", rownames(rna_sig_d1)),]
  rna_mat_d1 <- rna_mat_d1[-grep("chrX", rownames(rna_mat_d1)),]
  rna_pat_d1 <- rna_pat_d1[-grep("chrX", rownames(rna_pat_d1)),]
}



#rna_mat_d1 <- rna_mat_d1[rna_sig_d1$pvals_summed<=0.05,]
#rna_pat_d1 <- rna_pat_d1[rna_sig_d1$pvals_summed<=0.05,]
#rna_sig_d1 <- rna_sig_d1[rna_sig_d1$pvals_summed<=0.05,]

rna_mat_d1 <- rna_mat_d1[rna_sig_d1$pvals_summed<=0.001,]
rna_pat_d1 <- rna_pat_d1[rna_sig_d1$pvals_summed<=0.001,]
rna_sig_d1 <- rna_sig_d1[rna_sig_d1$pvals_summed<=0.001,]
rna_sig_d1_regions <- matrix(nrow=nrow(rna_sig_d1), ncol=3, data=unlist(strsplit(rownames(rna_sig_d1), split=":|-")), byrow=T)

################################################################################
#Read in phase groupings and identify the phase of each variant.
################################################################################

phase_grouping_d1 <- read.table(phase_grouping_d1, header=F, sep="\t", stringsAsFactors=F)


rna_sig_d1_regions_wPhase <- getPhaseGroup(rna_sig_d1_regions, phase_grouping_d1)




################################################################################
#Read in the TF tables. Remove X chrom. Remove histones, POL.
#Restrict to significant.
################################################################################


tf_sig_d1 <- read.table(tf_sig_d1, header=T, sep="\t", stringsAsFactors=F)


tf_mat_d1 <- read.table(tf_mat_d1, header=T, sep="\t", stringsAsFactors=F)

tf_pat_d1 <- read.table(tf_pat_d1, header=T, sep="\t", stringsAsFactors=F)


forbiddenExprs <- colnames(tf_sig_d1)[c(grep("H[3|4]K", colnames(tf_sig_d1)), grep("INPUT", colnames(tf_sig_d1)), grep("^POL", colnames(tf_sig_d1)))]
tf_sig_d1 <- tf_sig_d1[,!colnames(tf_sig_d1)%in%forbiddenExprs]
tf_mat_d1 <- tf_mat_d1[,!colnames(tf_mat_d1)%in%forbiddenExprs]
tf_pat_d1 <- tf_pat_d1[,!colnames(tf_pat_d1)%in%forbiddenExprs]

toRemove <- grep("chrX", rownames(tf_sig_d1))
if(length(toRemove)>0) {
  tf_sig_d1 <- tf_sig_d1[-grep("chrX", rownames(tf_sig_d1)),]
  tf_mat_d1 <- tf_mat_d1[-grep("chrX", rownames(tf_mat_d1)),]
  tf_pat_d1 <- tf_pat_d1[-grep("chrX", rownames(tf_pat_d1)),]

}


tf_sig_d1_minP <- rowMins(as.matrix(tf_sig_d1))
#tf_sig_d1 <- tf_sig_d1[tf_sig_d1_minP<=0.05,]
tf_sig_d1 <- tf_sig_d1[tf_sig_d1_minP<=0.001,]
tf_mat_d1 <- tf_mat_d1[tf_sig_d1_minP<=0.001,]
tf_pat_d1 <- tf_pat_d1[tf_sig_d1_minP<=0.001,]

tf_sig_d1_regions <- matrix(nrow=nrow(tf_sig_d1), ncol=3, data=unlist(strsplit(rownames(tf_sig_d1), split=":|-")), byrow=T)

################################################################################
#Determine phase of TF variants.
################################################################################

tf_sig_d1_regions_wPhase <- getPhaseGroup(tf_sig_d1_regions, phase_grouping_d1)


################################################################################
#Identify the cases of RNA and TF both in phase..
################################################################################


saveFile <- paste(outDir, "rna_sig_d1_regions_wPhase_w_SignificantTFInPhase.txt", sep="")

rna_sig_d1_regions_wPhase_w_SignificantTFInPhase <- getPhasedTFCases(rna_sig_d1_regions_wPhase, tf_sig_d1_regions_wPhase, saveFile)
#rna_sig_d1_regions_wPhase_w_SignificantTFInPhase <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

################################################################################
#For each case, identify the RNA pvalue.
#Include two columns which give the RNA mat and pat reads.
#For each case, identify all TFs affected.
#Include a column which, for all affected TFs, gives the pvalue (in order)
#Include two columns which, for all affected TFs, gives the reads for mat and apt.
################################################################################

infoTable <- c()
for (i in 1:nrow(rna_sig_d1_regions_wPhase_w_SignificantTFInPhase)) {
  if(i%%100==0) {print(paste(i, nrow(rna_sig_d1_regions_wPhase_w_SignificantTFInPhase)))}
  thisRNA_region <- paste(rna_sig_d1_regions_wPhase_w_SignificantTFInPhase[i,1], ":", rna_sig_d1_regions_wPhase_w_SignificantTFInPhase[i,2], "-", rna_sig_d1_regions_wPhase_w_SignificantTFInPhase[i,3], sep="")
  thisTF_region <- paste(rna_sig_d1_regions_wPhase_w_SignificantTFInPhase[i,4], ":", rna_sig_d1_regions_wPhase_w_SignificantTFInPhase[i,5], "-", rna_sig_d1_regions_wPhase_w_SignificantTFInPhase[i,6], sep="")

  rna_pval <- rna_sig_d1[thisRNA_region,"pvals_summed"]
  rna_matReads <- rna_mat_d1[thisRNA_region,"matSums"]
  rna_patReads <- rna_pat_d1[thisRNA_region,"patSums"]

  tf_names <- colnames(tf_sig_d1)[which(as.numeric(tf_sig_d1[thisTF_region,])<=0.001)]
  tf_pval <- tf_sig_d1[thisTF_region,tf_names]
  tf_matReads <- tf_mat_d1[thisTF_region,tf_names]
  tf_patReads <- tf_pat_d1[thisTF_region,tf_names]
  numTFs <- length(tf_names)
  if(length(tf_names)>1) {
    tf_names <- paste(tf_names, collapse=";")
    tf_pval <- paste(tf_pval, collapse=";")
    tf_matReads <- paste(tf_matReads, collapse=";")
    tf_patReads <- paste(tf_patReads, collapse=";")
  }
  thisLine <- c(rna_pval, rna_matReads, rna_patReads, numTFs, tf_names, tf_pval, tf_matReads, tf_patReads)
  infoTable <- rbind(infoTable, thisLine)
}

infoTable <- as.data.frame(infoTable)
rownames(infoTable) <- c()
colnames(infoTable) <- c("rna_pval", "rna_matReads", "rna_patReads", "numTFs", "tf_names", "tf_pval", "tf_matReads", "tf_patReads")
for (i in 1:4) { infoTable[,i] <- as.numeric(as.character(infoTable[,i]))}


rna_sig_d1_regions_wPhase_w_SignificantTFInPhase <- cbind(rna_sig_d1_regions_wPhase_w_SignificantTFInPhase, infoTable)

saveFile <- paste(outDir, "rna_sig_d1_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
write.table(rna_sig_d1_regions_wPhase_w_SignificantTFInPhase, saveFile, row.names=F, col.names=T, sep="\t", quote=F)




################################################################################
################################################################################
#Donor 2
################################################################################
################################################################################



################################################################################
#Read in the RNA tables. Remove X chrom. Restrict to significant.
################################################################################


rna_sig_d2 <- read.table(rna_sig_d2, header=T, sep="\t", stringsAsFactors=F)
rna_mat_d2 <- read.table(rna_mat_d2, header=T, sep="\t", stringsAsFactors=F)
rna_pat_d2 <- read.table(rna_pat_d2, header=T, sep="\t", stringsAsFactors=F)


toRemove <- grep("chrX", rownames(rna_sig_d2))
if(length(toRemove)>0) {
  rna_sig_d2 <- rna_sig_d2[-grep("chrX", rownames(rna_sig_d2)),]
  rna_mat_d2 <- rna_mat_d2[-grep("chrX", rownames(rna_mat_d2)),]
  rna_pat_d2 <- rna_pat_d2[-grep("chrX", rownames(rna_pat_d2)),]
}


#rna_mat_d2 <- rna_mat_d2[rna_sig_d2$pvals_summed<=0.05,]
#rna_pat_d2 <- rna_pat_d2[rna_sig_d2$pvals_summed<=0.05,]
#rna_sig_d2 <- rna_sig_d2[rna_sig_d2$pvals_summed<=0.05,]

rna_mat_d2 <- rna_mat_d2[rna_sig_d2$pvals_summed<=0.001,]
rna_pat_d2 <- rna_pat_d2[rna_sig_d2$pvals_summed<=0.001,]
rna_sig_d2 <- rna_sig_d2[rna_sig_d2$pvals_summed<=0.001,]
rna_sig_d2_regions <- matrix(nrow=nrow(rna_sig_d2), ncol=3, data=unlist(strsplit(rownames(rna_sig_d2), split=":|-")), byrow=T)



################################################################################
#Read in phase groupings and identify the phase of each variant.
################################################################################


phase_grouping_d2 <- read.table(phase_grouping_d2, header=F, sep="\t", stringsAsFactors=F)


rna_sig_d2_regions_wPhase <- getPhaseGroup(rna_sig_d2_regions, phase_grouping_d2)


################################################################################
#Read in the TF tables. Remove X chrom. Remove histones, POL.
#Restrict to significant.
################################################################################


tf_sig_d2 <- read.table(tf_sig_d2, header=T, sep="\t", stringsAsFactors=F)

tf_mat_d2 <- read.table(tf_mat_d2, header=T, sep="\t", stringsAsFactors=F)

tf_pat_d2 <- read.table(tf_pat_d2, header=T, sep="\t", stringsAsFactors=F)



forbiddenExprs <- colnames(tf_sig_d2)[c(grep("H[3|4]K", colnames(tf_sig_d2)), grep("INPUT", colnames(tf_sig_d2)), grep("^POL", colnames(tf_sig_d2)))]
tf_sig_d2 <- tf_sig_d2[,!colnames(tf_sig_d2)%in%forbiddenExprs]
tf_mat_d2 <- tf_mat_d2[,!colnames(tf_mat_d2)%in%forbiddenExprs]
tf_pat_d2 <- tf_pat_d2[,!colnames(tf_pat_d2)%in%forbiddenExprs]

toRemove <- grep("chrX", rownames(tf_sig_d2))
if(length(toRemove)>0) {
  tf_sig_d2 <- tf_sig_d2[-grep("chrX", rownames(tf_sig_d2)),]
  tf_mat_d2 <- tf_mat_d2[-grep("chrX", rownames(tf_mat_d2)),]
  tf_pat_d2 <- tf_pat_d2[-grep("chrX", rownames(tf_pat_d2)),]
}


tf_sig_d2_minP <- rowMins(as.matrix(tf_sig_d2))
tf_sig_d2 <- tf_sig_d2[tf_sig_d2_minP<=0.001,]
tf_mat_d2 <- tf_mat_d2[tf_sig_d2_minP<=0.001,]
tf_pat_d2 <- tf_pat_d2[tf_sig_d2_minP<=0.001,]

tf_sig_d2_regions <- matrix(nrow=nrow(tf_sig_d2), ncol=3, data=unlist(strsplit(rownames(tf_sig_d2), split=":|-")), byrow=T)



################################################################################
#Determine phase of TF variants.
################################################################################

tf_sig_d2_regions_wPhase <- getPhaseGroup(tf_sig_d2_regions, phase_grouping_d2)


################################################################################
#Identify the cases of RNA and TF both in phase..
################################################################################


saveFile <- paste(outDir, "rna_sig_d2_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
rna_sig_d2_regions_wPhase_w_SignificantTFInPhase <- getPhasedTFCases(rna_sig_d2_regions_wPhase, tf_sig_d2_regions_wPhase, saveFile)
#rna_sig_d2_regions_wPhase_w_SignificantTFInPhase <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)



################################################################################
#For each case, identify the RNA pvalue.
#Include two columns which give the RNA mat and pat reads.
#For each case, identify all TFs affected.
#Include a column which, for all affected TFs, gives the pvalue (in order)
#Include two columns which, for all affected TFs, gives the reads for mat and apt.
################################################################################

infoTable <- c()
for (i in 1:nrow(rna_sig_d2_regions_wPhase_w_SignificantTFInPhase)) {
  if(i%%100==0) {print(paste(i, nrow(rna_sig_d2_regions_wPhase_w_SignificantTFInPhase)))}
  thisRNA_region <- paste(rna_sig_d2_regions_wPhase_w_SignificantTFInPhase[i,1], ":", rna_sig_d2_regions_wPhase_w_SignificantTFInPhase[i,2], "-", rna_sig_d2_regions_wPhase_w_SignificantTFInPhase[i,3], sep="")
  thisTF_region <- paste(rna_sig_d2_regions_wPhase_w_SignificantTFInPhase[i,4], ":", rna_sig_d2_regions_wPhase_w_SignificantTFInPhase[i,5], "-", rna_sig_d2_regions_wPhase_w_SignificantTFInPhase[i,6], sep="")

  rna_pval <- rna_sig_d2[thisRNA_region,"pvals_summed"]
  rna_matReads <- rna_mat_d2[thisRNA_region,"matSums"]
  rna_patReads <- rna_pat_d2[thisRNA_region,"patSums"]

  tf_names <- colnames(tf_sig_d2)[which(as.numeric(tf_sig_d2[thisTF_region,])<=0.001)]
  tf_pval <- tf_sig_d2[thisTF_region,tf_names]
  tf_matReads <- tf_mat_d2[thisTF_region,tf_names]
  tf_patReads <- tf_pat_d2[thisTF_region,tf_names]
  numTFs <- length(tf_names)
  if(length(tf_names)>1) {
    tf_names <- paste(tf_names, collapse=";")
    tf_pval <- paste(tf_pval, collapse=";")
    tf_matReads <- paste(tf_matReads, collapse=";")
    tf_patReads <- paste(tf_patReads, collapse=";")
  }
  thisLine <- c(rna_pval, rna_matReads, rna_patReads, numTFs, tf_names, tf_pval, tf_matReads, tf_patReads)
  infoTable <- rbind(infoTable, thisLine)
}

infoTable <- as.data.frame(infoTable)
rownames(infoTable) <- c()
colnames(infoTable) <- c("rna_pval", "rna_matReads", "rna_patReads", "numTFs", "tf_names", "tf_pval", "tf_matReads", "tf_patReads")
for (i in 1:4) { infoTable[,i] <- as.numeric(as.character(infoTable[,i]))}


rna_sig_d2_regions_wPhase_w_SignificantTFInPhase <- cbind(rna_sig_d2_regions_wPhase_w_SignificantTFInPhase, infoTable)

saveFile <- paste(outDir, "rna_sig_d2_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
write.table(rna_sig_d2_regions_wPhase_w_SignificantTFInPhase, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
