#!/usr/bin/env R
#Figure_2A_2B.R



################################################################################
#This script is used to produce Figures 2A, 2B, and supplemental figures .
#
#Run under R 4.1.0, but can be used under other versions.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# varSeqs_d1: A FASTA file for donor 1 which contains all regions that contain a
#     non-reference base, whether heterozygous or homozygous. Produced by the
#     Building_fasta_variant_sequences.sh script.
# varSeqs_d2: as varSeqs_d1, but for donor2.
# pval_table_1: A table compiling p-value calculations across all experiments for
#     donor 1, produced by the script CompilingHaplotypePvalsTable.R
# pval_table_2: As pval_table_1, but for donor 2.
# depthTable_1_mat: A table compiling the number of reads corresponding to the
#     number of reads for each variant region in each dataset for haplotype 1.
#     Produced by the script CompilingHaplotypeReadsTable.R .
# depthTable_1_pat: As depthTable_1_mat, but for the alternate haplotype.
# depthTable_2_mat: As depthTable_1_mat, but for donor 2.
# depthTable_2_pat: As depthTable_2_mat, but for the alternate haplotype.
#
################################################################################



################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################



library(GenomicRanges)
library(ggplot2)
library(matrixStats)


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
varSeqs_d1 <- args[2]
varSeqs_d2 <- args[3]
pval_table_1 <- args[4]
pval_table_2 <- args[5]
depthTable_1_mat <- args[6]
depthTable_1_pat <- args[7]
depthTable_2_mat <- args[8]
depthTable_2_pat <- args[9]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#varSeqs_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#varSeqs_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg.txt"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg.txt"
#depthTable_1_mat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg.txt"
#depthTable_1_pat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg.txt"
#depthTable_2_mat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg.txt"
#depthTable_2_pat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg.txt"




################################################################################
#Load in the shared regions across donors first.
#Make sure we remove files that are no longer necessary once analysis is complete.
################################################################################

varSeqs_d1 <- readLines(varSeqs_d1)
varSeqs_d1 <- varSeqs_d1[grep("^>", varSeqs_d1)]
varSeqs_d1_m <- varSeqs_d1[grep(":M:", varSeqs_d1)]
varSeqs_d1_m <- matrix(nrow=length(varSeqs_d1_m), ncol=2, data=unlist(strsplit(varSeqs_d1_m, split=":M:")), byrow=T)
varSeqs_d1_p <- varSeqs_d1[grep(":P:", varSeqs_d1)]
varSeqs_d1_p <- matrix(nrow=length(varSeqs_d1_p), ncol=2, data=unlist(strsplit(varSeqs_d1_p, split=":P:")), byrow=T)
varSeqs_d1 <- cbind(varSeqs_d1_m, varSeqs_d1_p[,2])
colnames(varSeqs_d1) <- c("Region", "HapM", "HapP")
varSeqs_d1[,"Region"] <- gsub("^>", "", varSeqs_d1[,"Region"])

varSeqs_d2 <- readLines(varSeqs_d2)
varSeqs_d2 <- varSeqs_d2[grep("^>", varSeqs_d2)]
varSeqs_d2_m <- varSeqs_d2[grep(":M:", varSeqs_d2)]
varSeqs_d2_m <- matrix(nrow=length(varSeqs_d2_m), ncol=2, data=unlist(strsplit(varSeqs_d2_m, split=":M:")), byrow=T)
varSeqs_d2_p <- varSeqs_d2[grep(":P:", varSeqs_d2)]
varSeqs_d2_p <- matrix(nrow=length(varSeqs_d2_p), ncol=2, data=unlist(strsplit(varSeqs_d2_p, split=":P:")), byrow=T)
varSeqs_d2 <- cbind(varSeqs_d2_m, varSeqs_d2_p[,2])
colnames(varSeqs_d2) <- c("Region", "HapM", "HapP")
varSeqs_d2[,"Region"] <- gsub("^>", "", varSeqs_d2[,"Region"])


toRemove <- grep("chrX", varSeqs_d1[,1])
if(length(toRemove)>0) {
  varSeqs_d1 <- varSeqs_d1[-toRemove,]
}

toRemove <- grep("chrX", varSeqs_d2[,1])
if(length(toRemove)>0) {
  varSeqs_d2 <- varSeqs_d2[-toRemove,]
}


varSeqs_d1_shared <- varSeqs_d1[varSeqs_d1[,"Region"]%in%varSeqs_d2[,"Region"],]
varSeqs_d2_shared <- varSeqs_d2[varSeqs_d2[,"Region"]%in%varSeqs_d1[,"Region"],]

varSeqs_d1_shared <- varSeqs_d1_shared[order(varSeqs_d1_shared[,"Region"]),]
varSeqs_d2_shared <- varSeqs_d2_shared[order(varSeqs_d2_shared[,"Region"]),]

shared_variant_matrix <- matrix(nrow=nrow(varSeqs_d1_shared), ncol=6, data=NA)
for (i in 1:nrow(varSeqs_d1_shared)) {
  if(varSeqs_d1_shared[i,1]==varSeqs_d2_shared[i,1] && varSeqs_d1_shared[i,2]==varSeqs_d2_shared[i,2]) {
    shared_variant_matrix[i,1] <- varSeqs_d1_shared[i,"Region"]
    shared_variant_matrix[i,2] <- 1
    shared_variant_matrix[i,3] <- varSeqs_d1_shared[i,"HapM"]
    shared_variant_matrix[i,4] <- varSeqs_d1_shared[i,"HapP"]
    shared_variant_matrix[i,5] <- varSeqs_d2_shared[i,"HapM"]
    shared_variant_matrix[i,6] <- varSeqs_d2_shared[i,"HapP"]
  }
  if(varSeqs_d1_shared[i,1]==varSeqs_d2_shared[i,2] && varSeqs_d1_shared[i,2]==varSeqs_d2_shared[i,1]) {
    shared_variant_matrix[i,1] <- varSeqs_d1_shared[i,"Region"]
    shared_variant_matrix[i,2] <- 2
    shared_variant_matrix[i,3] <- varSeqs_d1_shared[i,"HapM"]
    shared_variant_matrix[i,4] <- varSeqs_d1_shared[i,"HapP"]
    shared_variant_matrix[i,5] <- varSeqs_d2_shared[i,"HapM"]
    shared_variant_matrix[i,6] <- varSeqs_d2_shared[i,"HapP"]
  }
}
shared_variant_matrix <- shared_variant_matrix[!is.na(shared_variant_matrix[,1]),]





################################################################################
#Next, load in pval tables and depth tables.  Restrict
#these tables to only relevant regions ASAP.
################################################################################

pval_table_1 <- read.table(pval_table_1, header=T, sep="\t", stringsAsFactors=F)

sigInput_d1 <- pval_table_1[which(rowMins(as.matrix(pval_table_1[grep("INPUT", colnames(pval_table_1))]))<=0.05),]
sigInput_d1 <- sigInput_d1[-grep("chrX", rownames(sigInput_d1)),]

allMat <- matrix(nrow=ncol(pval_table_1), ncol=3, data=unlist(strsplit(colnames(pval_table_1), split="_")), byrow=T)
allTFs <- unique(allMat[,2])
allTissues <- unique(allMat[,1])

forbiddenExprs <- c("CB_ZNF207_1224", "CB_ASCL1_1224", "CB_ATF7_1224", paste(allTissues, "OLIG2_1224", sep="_"), "CB_OLIG2_1230", "CB_ASCL1_1224")
if(length(forbiddenExprs[forbiddenExprs%in%colnames(pval_table_1)])>0) { pval_table_1 <- pval_table_1[,-which(colnames(pval_table_1)%in%forbiddenExprs)]}

forbiddenExprs_2 <- c(colnames(pval_table_1)[grep("H3K", colnames(pval_table_1))], colnames(pval_table_1)[grep("INPUT", colnames(pval_table_1))], colnames(pval_table_1)[grep("POL2", colnames(pval_table_1))])
pval_table_1 <- pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs_2]

depthTable_1_mat <- read.table(depthTable_1_mat, header=T, sep="\t", stringsAsFactors=F)
depthTable_1_pat <- read.table(depthTable_1_pat, header=T, sep="\t", stringsAsFactors=F)

depthTable_1_mat <- depthTable_1_mat[,colnames(depthTable_1_mat)%in%colnames(pval_table_1)]
depthTable_1_pat <- depthTable_1_pat[,colnames(depthTable_1_pat)%in%colnames(pval_table_1)]

toRemove <- grep("chrX", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
  depthTable_1_mat <- depthTable_1_mat[-toRemove,]
  depthTable_1_pat <- depthTable_1_pat[-toRemove,]
}

shared_pval_table_1 <- pval_table_1[rownames(pval_table_1)%in%shared_variant_matrix[,1],]
shared_depthTable_1_mat <- depthTable_1_mat[rownames(depthTable_1_mat)%in%shared_variant_matrix[,1],]
shared_depthTable_1_pat <- depthTable_1_pat[rownames(depthTable_1_pat)%in%shared_variant_matrix[,1],]







pval_table_2 <- read.table(pval_table_2, header=T, sep="\t", stringsAsFactors=F)

sigInput_d2 <- pval_table_2[which(rowMins(as.matrix(pval_table_2[grep("INPUT", colnames(pval_table_2))]))<=0.05),]

forbiddenExprs <- c("CB_ZNF207_1224", "CB_ASCL1_1224", "CB_ATF7_1224", paste(allTissues, "OLIG2_1224", sep="_"), "CB_OLIG2_1230", "CB_ASCL1_1224")
if(length(forbiddenExprs[forbiddenExprs%in%colnames(pval_table_2)])>0) { pval_table_2 <- pval_table_2[,-which(colnames(pval_table_2)%in%forbiddenExprs)]}

forbiddenExprs_2 <- c(colnames(pval_table_2)[grep("H3K", colnames(pval_table_2))], colnames(pval_table_2)[grep("INPUT", colnames(pval_table_2))], colnames(pval_table_2)[grep("POL2", colnames(pval_table_2))])
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs_2]

depthTable_2_mat <- read.table(depthTable_2_mat, header=T, sep="\t", stringsAsFactors=F)
depthTable_2_pat <- read.table(depthTable_2_pat, header=T, sep="\t", stringsAsFactors=F)

depthTable_2_mat <- depthTable_2_mat[,colnames(depthTable_2_mat)%in%colnames(pval_table_2)]
depthTable_2_pat <- depthTable_2_pat[,colnames(depthTable_2_pat)%in%colnames(pval_table_2)]

toRemove <- grep("chrX", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
  depthTable_2_mat <- depthTable_2_mat[-toRemove,]
  depthTable_2_pat <- depthTable_2_pat[-toRemove,]
}

allMat <- matrix(nrow=ncol(pval_table_2), ncol=3, data=unlist(strsplit(colnames(pval_table_2), split="_")), byrow=T)
allTFs <- unique(allMat[,2])
allTissues <- unique(allMat[,1])



forbiddenExprs <- c("CB_ZNF207_1224", "CB_ASCL1_1224", "CB_ATF7_1224", paste(allTissues, "OLIG2_1224", sep="_"), "CB_OLIG2_1230", "CB_ASCL1_1224")

pval_table_1 <- pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]
depthTable_1_mat <- depthTable_1_mat[,!colnames(depthTable_1_mat)%in%forbiddenExprs]
depthTable_1_pat <- depthTable_1_pat[,!colnames(depthTable_1_pat)%in%forbiddenExprs]
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs]
depthTable_2_mat <- depthTable_2_mat[,!colnames(depthTable_2_mat)%in%forbiddenExprs]
depthTable_2_pat <- depthTable_2_pat[,!colnames(depthTable_2_pat)%in%forbiddenExprs]





colnames(pval_table_1) <- gsub("_1224", "", colnames(pval_table_1))
colnames(depthTable_1_mat) <- gsub("_1224", "", colnames(depthTable_1_mat))
colnames(depthTable_1_pat) <- gsub("_1224", "", colnames(depthTable_1_pat))

colnames(pval_table_2) <- gsub("_1230", "", colnames(pval_table_2))
colnames(depthTable_2_mat) <- gsub("_1230", "", colnames(depthTable_2_mat))
colnames(depthTable_2_pat) <- gsub("_1230", "", colnames(depthTable_2_pat))







shared_pval_table_1 <- pval_table_1[,colnames(pval_table_1)%in%colnames(pval_table_2)]
shared_depthTable_1_mat <- depthTable_1_mat[,colnames(depthTable_1_mat)%in%colnames(depthTable_2_mat)]
shared_depthTable_1_pat <- depthTable_1_pat[,colnames(depthTable_1_pat)%in%colnames(depthTable_2_pat)]

shared_pval_table_2 <- pval_table_2[,colnames(pval_table_2)%in%colnames(pval_table_1)]
shared_depthTable_2_mat <- depthTable_2_mat[,colnames(depthTable_2_mat)%in%colnames(depthTable_1_mat)]
shared_depthTable_2_pat <- depthTable_2_pat[,colnames(depthTable_2_pat)%in%colnames(depthTable_1_pat)]

shared_pval_table_1 <- shared_pval_table_1[rownames(shared_pval_table_1)%in%shared_variant_matrix[,1],]
shared_depthTable_1_mat <- shared_depthTable_1_mat[rownames(shared_depthTable_1_mat)%in%shared_variant_matrix[,1],]
shared_depthTable_1_pat <- shared_depthTable_1_pat[rownames(shared_depthTable_1_pat)%in%shared_variant_matrix[,1],]

shared_pval_table_2 <- shared_pval_table_2[rownames(shared_pval_table_2)%in%shared_variant_matrix[,1],]
shared_depthTable_2_mat <- shared_depthTable_2_mat[rownames(shared_depthTable_2_mat)%in%shared_variant_matrix[,1],]
shared_depthTable_2_pat <- shared_depthTable_2_pat[rownames(shared_depthTable_2_pat)%in%shared_variant_matrix[,1],]

shared_pval_table_1 <- shared_pval_table_1[,order(colnames(shared_pval_table_1))]
shared_depthTable_1_mat <- shared_depthTable_1_mat[,order(colnames(shared_depthTable_1_mat))]
shared_depthTable_1_pat <- shared_depthTable_1_pat[,order(colnames(shared_depthTable_1_pat))]
shared_pval_table_2 <- shared_pval_table_2[,order(colnames(shared_pval_table_2))]
shared_depthTable_2_mat <- shared_depthTable_2_mat[,order(colnames(shared_depthTable_2_mat))]
shared_depthTable_2_pat <- shared_depthTable_2_pat[,order(colnames(shared_depthTable_2_pat))]

allMat_shared <- matrix(nrow=ncol(shared_pval_table_1), ncol=2, data=unlist(strsplit(colnames(shared_pval_table_1), split="_")), byrow=T)
allTFs_shared <- unique(allMat_shared[,2])
allTissues_shared <- unique(allMat_shared[,1])








################################################################################
#Set up the tables with the noInput cases.
################################################################################

shared_variant_matrix_noInput <- shared_variant_matrix[!shared_variant_matrix[,1]%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]

shared_pval_table_1_noInput <- shared_pval_table_1[!rownames(shared_pval_table_1)%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]
shared_depthTable_1_mat_noInput <- shared_depthTable_1_mat[!rownames(shared_depthTable_1_mat)%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]
shared_depthTable_1_pat_noInput <- shared_depthTable_1_pat[!rownames(shared_depthTable_1_pat)%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]
shared_pval_table_2_noInput <- shared_pval_table_2[!rownames(shared_pval_table_2)%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]
shared_depthTable_2_mat_noInput <- shared_depthTable_2_mat[!rownames(shared_depthTable_2_mat)%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]
shared_depthTable_2_pat_noInput <- shared_depthTable_2_pat[!rownames(shared_depthTable_2_pat)%in%c(rownames(sigInput_d1), rownames(sigInput_d2)),]


################################################################################
#First, across-donor consistency.
#For each TF in each tissue, for each donor, identify significant variants.
#Given all significant variants, determine whether or not the bias is in the
#same direction-- keeping in mind that there may be a PAT/MAT assignment switch
#(Designated as 1 for a match and 2 for a mismatch)
################################################################################

################################################################################
#Identify p-value cutoffs and experiments
################################################################################

forbiddenFactors <- c(allTFs[grep("H3K", allTFs)], allTFs[grep("INPUT", allTFs)], allTFs[grep("POL2", allTFs)])
allTFs <- allTFs[!allTFs%in%forbiddenFactors]

pval_cutoffs <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001,
  0.000000001, 0.0000000001)




################################################################################
#Run the checks.
################################################################################

donor_consistency_table_noInput <- c()
for (i in 1:ncol(shared_pval_table_1_noInput)) {
  print(paste(i, ncol(shared_pval_table_1_noInput)))
  thisPval_d1 <- shared_pval_table_1_noInput[,i]
  thisPval_d2 <- shared_pval_table_2_noInput[,1]

  thisMat_d1 <- shared_depthTable_1_mat_noInput[,1]
  thisMat_d2 <- shared_depthTable_2_mat_noInput[,1]

  thisPat_d1 <- shared_depthTable_1_pat_noInput[,1]
  thisPat_d2 <- shared_depthTable_2_pat_noInput[,1]

  this_totalDepth_d1 <- thisMat_d1 + thisPat_d1
  this_totalDepth_d2 <- thisMat_d2 + thisPat_d2

  theSigs_1 <- which(thisPval_d1<=0.05)
  theSigs_2 <- which(thisPval_d2<=0.05)
  theSigs_all <- unique(c(theSigs_1, theSigs_2))
  theSigs_all <- theSigs_all[order(theSigs_all)]
  minPvals <- rowMins(as.matrix(cbind(thisPval_d1[theSigs_all], thisPval_d2[theSigs_all])))

  this_consistency_table <- c()

  for (j in 1:length(theSigs_all)) {
    match_orientation <- as.numeric(shared_variant_matrix_noInput[theSigs_all[j],2])

    direction_1 <- "EQUAL"
    if(thisMat_d1[theSigs_all[j]]>thisPat_d1[theSigs_all[j]]) { direction_1 <- "MAT"}
    if(thisMat_d1[theSigs_all[j]]<thisPat_d1[theSigs_all[j]]) { direction_1 <- "PAT"}


    direction_2 <- "EQUAL"
    if(thisMat_d2[theSigs_all[j]]>thisPat_d2[theSigs_all[j]]) { direction_2 <- "MAT"}
    if(thisMat_d2[theSigs_all[j]]<thisPat_d2[theSigs_all[j]]) { direction_2 <- "PAT"}

    if(this_totalDepth_d1[theSigs_all[j]]>=6 && this_totalDepth_d2[theSigs_all[j]]>=6) {
      consistency <- "Inconsistent"
      if(match_orientation==1) {
        if(direction_1=="MAT" && direction_2=="MAT") { consistency <- "Consistent"}
        if(direction_1=="PAT" && direction_2=="PAT") { consistency <- "Consistent"}
      }
      if(match_orientation==2) {
        if(direction_1=="PAT" && direction_2=="MAT") { consistency <- "Consistent"}
        if(direction_1=="MAT" && direction_2=="PAT") { consistency <- "Consistent"}
      }

      thisVariant <- shared_variant_matrix_noInput[theSigs_all[j],1]
      thisLine <- c(thisVariant, allTFs[i], minPvals[j], direction_1, direction_2, consistency)
      this_consistency_table <- rbind(this_consistency_table, thisLine)
    }
  }
  this_consistency_table <- as.data.frame(this_consistency_table)
  colnames(this_consistency_table) <- c("variant", "tf", "pval", "direction_1", "direction_2", "consistency")
  this_consistency_table$pval <- as.numeric(as.character(this_consistency_table$pval))

  for (j in 1:length(pval_cutoffs)){
    thisSet <- this_consistency_table[this_consistency_table$pval<=pval_cutoffs[j],]
    numConsistent <- nrow(thisSet[thisSet$consistency=="Consistent",])
    numInconsistent <- nrow(thisSet[thisSet$consistency=="Inconsistent",])
    fractionConsistent <- numConsistent/nrow(thisSet)
    fractionInconsistent <- numInconsistent/nrow(thisSet)
    thisLine <- c(allTFs[i], pval_cutoffs[j], numConsistent, numInconsistent, fractionConsistent, fractionInconsistent)
    donor_consistency_table_noInput <- rbind(donor_consistency_table_noInput, thisLine)
  }

}


donor_consistency_table_noInput <- as.data.frame(donor_consistency_table_noInput)
colnames(donor_consistency_table_noInput) <- c("TF", "pvalue_cutoff", "numConsistent", "numInconsistent", "fractionConsistent", "fractionInconsistent")
for (i in 2:ncol(donor_consistency_table_noInput)) {donor_consistency_table_noInput[,i] <- as.numeric(as.character(donor_consistency_table_noInput[,i]))}


#saveFile <- paste(outDir, "Consistency_of_variants_acrossDonors_directionOnly_through_pval_cutoffs_noInput.txt", sep="")
#write.table(donor_consistency_table_noInput, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#donor_consistency_table_noInput <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)





################################################################################
#Make a reproducibility figure.
################################################################################

donor_consistency_table_noInput$fractionReproducible <- 1 - (donor_consistency_table_noInput$fractionInconsistent)*2


graphDF <- c()
for (i in 1:length(pval_cutoffs)) {
  thisSet <- donor_consistency_table_noInput[donor_consistency_table_noInput$pvalue_cutoff==pval_cutoffs[i],]
  total_consistent <- sum(thisSet$numConsistent)
  total_inconsistent <- sum(thisSet$numInconsistent)
  fraction_consistent <- total_consistent/(total_consistent+total_inconsistent)
  fraction_inconsistent <- total_inconsistent/(total_consistent+total_inconsistent)
  fraction_reproducible <- 1 - (fraction_inconsistent*2)
  graphDF <- rbind(graphDF, c(pval_cutoffs[i], "Reproducible", fraction_reproducible))
}


colnames(graphDF) <- c("Pval_cutoff", "Reproducible", "Fraction")
graphDF <- as.data.frame(graphDF)

graphDF$Pval_cutoff <- factor(graphDF$Pval_cutoff, levels=pval_cutoffs)
graphDF$Fraction <- as.numeric(as.character(graphDF$Fraction))


labels <- paste(round(graphDF$Fraction*100, digits=2), "%", sep="")
graphDF$position <- graphDF$Fraction - 0.05


saveFile <- paste(outDir, "Figure_2A.pdf", sep="")
ggplot(graphDF, aes(y=Fraction, x=Pval_cutoff)) + geom_bar(stat="identity", fill="grey") + theme_classic(base_size=25) + theme(legend.position="null") +
  xlab("Pvalue cutoff") + ylab("Fraction") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = labels, y=position-0.075), size=8, angle=90) + ylim(0,1) + xlab("p-value cutoff") + ylab("Fraction Reproducible")
ggsave(saveFile)


################################################################################
#Now, across-tissue consistency.
#For each TF in a given donor, identify all the experiments with that TF.
#Next, for each experiment, identify the significant variants.
#For each variant, determine the number of paired experiments whose reads support
#the same direction of effect as the significant variant of interest.
#Once that's done, move on to the next experiment, and only consider novel pairings--
#not those significant in previous experiments.
################################################################################

forbiddenFactors <- c(allTFs[grep("H3K", allTFs)], allTFs[grep("INPUT", allTFs)], allTFs[grep("POL2", allTFs)])
allTFs <- allTFs[!allTFs%in%forbiddenFactors]

pval_cutoffs <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001,
  0.000000001, 0.0000000001)




################################################################################
#Check Within-donor consistency
################################################################################

forbiddenFactors <- c(allTFs[grep("H3K", allTFs)], allTFs[grep("INPUT", allTFs)], allTFs[grep("POL2", allTFs)])
allTFs <- allTFs[!allTFs%in%forbiddenFactors]

pval_cutoffs <- c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001,
  0.000000001, 0.0000000001)


pval_table_1_noInput <- pval_table_1[!rownames(pval_table_1)%in%rownames(sigInput_d1),]
depthTable_1_mat_noInput <- depthTable_1_mat[!rownames(depthTable_1_mat)%in%rownames(sigInput_d1),]
depthTable_1_pat_noInput <- depthTable_1_pat[!rownames(depthTable_1_pat)%in%rownames(sigInput_d1),]

tissue_consistency_d1_noInput <- c()
for (i in 1:length(allTFs)) {
  this_pvalSet <- pval_table_1_noInput[,grep(allTFs[i], colnames(pval_table_1_noInput))]
  this_matSet <- depthTable_1_mat_noInput[,grep(allTFs[i], colnames(depthTable_1_mat_noInput))]
  this_patSet <- depthTable_1_pat_noInput[,grep(allTFs[i], colnames(depthTable_1_pat_noInput))]

  this_pvalSet <- this_pvalSet[,order(colnames(this_pvalSet))]
  this_matSet <- this_matSet[,order(colnames(this_matSet))]
  this_patSet <- this_patSet[,order(colnames(this_patSet))]

  this_consistency_table <- c()

  for (j in 1:ncol(this_pvalSet)) {
    print(paste(i, length(allTFs), j, ncol(this_pvalSet), nrow(this_consistency_table)))
    theSigs_1 <- which(this_pvalSet[,j]<=0.05)
    tissue1 <- gsub(paste("_", allTFs[i], sep=""), "", colnames(this_pvalSet)[j])
    for (k in 1:ncol(this_pvalSet)){
      tissue2 <- gsub(paste("_", allTFs[i], sep=""), "", colnames(this_pvalSet)[k])
      if(j!=k) {
        theConsidered <- which(this_matSet[,k]+this_patSet[,k]>=6)
        theSigs_1 <- theSigs_1[theSigs_1%in%theConsidered]
        if(k<j) {
          theSigs_2 <- which(this_pvalSet[,k]<=0.05)
          theSigs_1 <- theSigs_1[!theSigs_1%in%theSigs_2]
        }
        for (l in 1:length(theSigs_1)) {
          thisVariant <- rownames(this_pvalSet)[theSigs_1[l]]
          pval1 <- this_pvalSet[theSigs_1[l],j]
          pval2 <- this_pvalSet[theSigs_1[l],k]
          direction_1 <- "EQUAL"
          if(this_matSet[theSigs_1[l],j]>this_patSet[theSigs_1[l],j]) { direction_1 <- "MAT"}
          if(this_matSet[theSigs_1[l],j]<this_patSet[theSigs_1[l],j]) { direction_1 <- "PAT"}
          direction_2 <- "EQUAL"
          if(this_matSet[theSigs_1[l],k]>this_patSet[theSigs_1[l],k]) { direction_2 <- "MAT"}
          if(this_matSet[theSigs_1[l],k]<this_patSet[theSigs_1[l],k]) { direction_2 <- "PAT"}
          consistency <- "Indeterminate"
          if(direction_1==direction_2) { consistency <- "Consistent"}
          if(direction_1!=direction_2) { consistency <- "Inconsistent"}

          thisLine <- c(thisVariant, allTFs[i], tissue1, tissue2, pval1, pval2, direction_1, direction_2, consistency)
          this_consistency_table <- rbind(this_consistency_table, thisLine)
        }
      }
    }
  }
  this_consistency_table <- as.data.frame(this_consistency_table)
  colnames(this_consistency_table) <- c("variant", "tf", "tissue1", "tissue2", "pval1", "pval2", "direction_1", "direction_2", "consistency")
  this_consistency_table$pval1 <- as.numeric(as.character(this_consistency_table$pval1))

  for (j in 1:length(pval_cutoffs)){
    thisSet <- this_consistency_table[this_consistency_table$pval1<=pval_cutoffs[j],]
    numConsistent <- nrow(thisSet[thisSet$consistency=="Consistent",])
    numInconsistent <- nrow(thisSet[thisSet$consistency=="Inconsistent",])
    fractionConsistent <- numConsistent/nrow(thisSet)
    fractionInconsistent <- numInconsistent/nrow(thisSet)
    thisLine <- c(allTFs[i], pval_cutoffs[j], numConsistent, numInconsistent, fractionConsistent, fractionInconsistent)
    tissue_consistency_d1_noInput <- rbind(tissue_consistency_d1_noInput, thisLine)
  }

}

tissue_consistency_d1_noInput <- as.data.frame(tissue_consistency_d1_noInput)
colnames(tissue_consistency_d1_noInput) <- c("TF", "pvalue_cutoff", "numConsistent", "numInconsistent", "fractionConsistent", "fractionInconsistent")
for (i in 2:ncol(tissue_consistency_d1_noInput)) {tissue_consistency_d1_noInput[,i] <- as.numeric(as.character(tissue_consistency_d1_noInput[,i]))}

#saveFile <- paste(outDir, "Consistency_of_variants_directionOnly_through_pval_cutoffs_d1_noInput.txt", sep="")
#write.table(tissue_consistency_d1_noInput, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#tissue_consistency_d1_noInput <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)





pval_table_2_noInput <- pval_table_2[!rownames(pval_table_2)%in%rownames(sigInput_d2),]
depthTable_2_mat_noInput <- depthTable_2_mat[!rownames(depthTable_2_mat)%in%rownames(sigInput_d2),]
depthTable_2_pat_noInput <- depthTable_2_pat[!rownames(depthTable_2_pat)%in%rownames(sigInput_d2),]

tissue_consistency_d2_noInput <- c()
for (i in 1:length(allTFs)) {
  this_pvalSet <- pval_table_2_noInput[,grep(allTFs[i], colnames(pval_table_2_noInput))]
  this_matSet <- depthTable_2_mat_noInput[,grep(allTFs[i], colnames(depthTable_2_mat_noInput))]
  this_patSet <- depthTable_2_pat_noInput[,grep(allTFs[i], colnames(depthTable_2_pat_noInput))]

  this_pvalSet <- this_pvalSet[,order(colnames(this_pvalSet))]
  this_matSet <- this_matSet[,order(colnames(this_matSet))]
  this_patSet <- this_patSet[,order(colnames(this_patSet))]

  this_consistency_table <- c()

  for (j in 1:ncol(this_pvalSet)) {
    print(paste(i, length(allTFs), j, ncol(this_pvalSet), nrow(this_consistency_table)))
    theSigs_1 <- which(this_pvalSet[,j]<=0.05)
    tissue1 <- gsub(paste("_", allTFs[i], sep=""), "", colnames(this_pvalSet)[j])
    for (k in 1:ncol(this_pvalSet)){
      tissue2 <- gsub(paste("_", allTFs[i], sep=""), "", colnames(this_pvalSet)[k])
      if(j!=k) {
        theConsidered <- which(this_matSet[,k]+this_patSet[,k]>=6)
        theSigs_1 <- theSigs_1[theSigs_1%in%theConsidered]
        if(k<j) {
          theSigs_2 <- which(this_pvalSet[,k]<=0.05)
          theSigs_1 <- theSigs_1[!theSigs_1%in%theSigs_2]
        }
        if(length(theSigs_1)>0) {
          for (l in 1:length(theSigs_1)) {
            thisVariant <- rownames(this_pvalSet)[theSigs_1[l]]
            pval1 <- this_pvalSet[theSigs_1[l],j]
            pval2 <- this_pvalSet[theSigs_1[l],k]
            direction_1 <- "EQUAL"
            if(this_matSet[theSigs_1[l],j]>this_patSet[theSigs_1[l],j]) { direction_1 <- "MAT"}
            if(this_matSet[theSigs_1[l],j]<this_patSet[theSigs_1[l],j]) { direction_1 <- "PAT"}
            direction_2 <- "EQUAL"
            if(this_matSet[theSigs_1[l],k]>this_patSet[theSigs_1[l],k]) { direction_2 <- "MAT"}
            if(this_matSet[theSigs_1[l],k]<this_patSet[theSigs_1[l],k]) { direction_2 <- "PAT"}
            consistency <- "Indeterminate"
            if(direction_1==direction_2) { consistency <- "Consistent"}
            if(direction_1!=direction_2) { consistency <- "Inconsistent"}

            thisLine <- c(thisVariant, allTFs[i], tissue1, tissue2, pval1, pval2, direction_1, direction_2, consistency)
            this_consistency_table <- rbind(this_consistency_table, thisLine)
          }
        }
      }
    }
  }
  this_consistency_table <- as.data.frame(this_consistency_table)
  colnames(this_consistency_table) <- c("variant", "tf", "tissue1", "tissue2", "pval1", "pval2", "direction_1", "direction_2", "consistency")
  this_consistency_table$pval1 <- as.numeric(as.character(this_consistency_table$pval1))

  for (j in 1:length(pval_cutoffs)){
    thisSet <- this_consistency_table[this_consistency_table$pval1<=pval_cutoffs[j],]
    numConsistent <- nrow(thisSet[thisSet$consistency=="Consistent",])
    numInconsistent <- nrow(thisSet[thisSet$consistency=="Inconsistent",])
    fractionConsistent <- numConsistent/nrow(thisSet)
    fractionInconsistent <- numInconsistent/nrow(thisSet)
    thisLine <- c(allTFs[i], pval_cutoffs[j], numConsistent, numInconsistent, fractionConsistent, fractionInconsistent)
    tissue_consistency_d2_noInput <- rbind(tissue_consistency_d2_noInput, thisLine)
  }

}

tissue_consistency_d2_noInput <- as.data.frame(tissue_consistency_d2_noInput)
colnames(tissue_consistency_d2_noInput) <- c("TF", "pvalue_cutoff", "numConsistent", "numInconsistent", "fractionConsistent", "fractionInconsistent")
for (i in 2:ncol(tissue_consistency_d2_noInput)) {tissue_consistency_d2_noInput[,i] <- as.numeric(as.character(tissue_consistency_d2_noInput[,i]))}


#saveFile <- paste(outDir, "Consistency_of_variants_directionOnly_through_pval_cutoffs_d2_noInput.txt", sep="")
#write.table(tissue_consistency_d2_noInput, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#tissue_consistency_d2_noInput <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)




graphDF <- c()
for (i in 1:length(pval_cutoffs)) {
  thisSet <- rbind(tissue_consistency_d1_noInput[tissue_consistency_d1_noInput$pvalue_cutoff==pval_cutoffs[i],], tissue_consistency_d2_noInput[tissue_consistency_d2_noInput$pvalue_cutoff==pval_cutoffs[i],])
  total_consistent <- sum(thisSet$numConsistent)
  total_inconsistent <- sum(thisSet$numInconsistent)
  fraction_consistent <- total_consistent/(total_consistent+total_inconsistent)
  fraction_inconsistent <- total_inconsistent/(total_consistent+total_inconsistent)
  fraction_reproducible <- 1 - (fraction_inconsistent*2)
  graphDF <- rbind(graphDF, c(pval_cutoffs[i], "Reproducible", fraction_reproducible))
}


colnames(graphDF) <- c("Pval_cutoff", "Reproducible", "Fraction")
graphDF <- as.data.frame(graphDF)

graphDF$Pval_cutoff <- factor(graphDF$Pval_cutoff, levels=pval_cutoffs)
graphDF$Fraction <- as.numeric(as.character(graphDF$Fraction))


labels <- paste(round(graphDF$Fraction*100, digits=2), "%", sep="")
graphDF$position <- graphDF$Fraction - 0.05


saveFile <- paste(outDir, "Figure_2B.pdf", sep="")
ggplot(graphDF, aes(y=Fraction, x=Pval_cutoff)) + geom_bar(stat="identity", fill="grey") + theme_classic(base_size=25) + theme(legend.position="null") +
  xlab("Pvalue cutoff") + ylab("Fraction") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = labels, y=position-0.075), size=8, angle=90) + ylim(0,1) + xlab("p-value cutoff") + ylab("Fraction Reproducible")
ggsave(saveFile)







#
