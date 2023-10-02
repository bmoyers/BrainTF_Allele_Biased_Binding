#!/usr/bin/env R
#Table_2.R



################################################################################
#This script is used to determine the numbers that went in to table 2.
#
#Run under R 4.1.0, but can be used under other versions so long as the packages
#     GenomicRanges, matrixStats, and ggplot are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# varSeqs_d1: path to the fasta file for all variant sequences in donor 1,
#     produced by the script Building_fasta_variant_sequences.sh
# varSeqs_d2: As varSeqs_d1, but for donor 2.
# pval_table_1_summed_saveFile: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# pval_table_2_summed_saveFile: As pval_table_1_summed_saveFile, but for donor 2
# pval_table_1_all_saveFile: A table compiling all p-values for each TF over all variants
#     for all datasets for donor 1, produced by the script CompilingHaplotypePvalsTable.R
# pval_table_2_all_saveFile: As pval_table_1_all_saveFile, but for donor 2.
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
pval_table_1_summed_saveFile <- args[4]
pval_table_2_summed_saveFile <- args[5]
pval_table_1_all_saveFile <- args[6]
pval_table_2_all_saveFile <- args[7]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/"
#varSeqs_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#varSeqs_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#pval_table_1_summed_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pval_table_2_summed_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pval_table_1_all_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg.txt"
#pval_table_2_all_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg.txt"




################################################################################
#Before Anything Else, load in the summed p-value tables and identify cases of
#significant input.
################################################################################

pval_table_1_summed <- read.delim(pval_table_1_summed_saveFile, header=T, sep="\t", stringsAsFactors=F)
pval_table_2_summed <- read.delim(pval_table_2_summed_saveFile, header=T, sep="\t", stringsAsFactors=F)

sigInput_1 <- rownames(pval_table_1_summed[which(pval_table_1_summed$INPUT<=0.05),])
sigInput_2 <- rownames(pval_table_2_summed[which(pval_table_2_summed$INPUT<=0.05),])


################################################################################
#Load in the shared regions and whatnot first.
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

dim(varSeqs_d1)
#3138109
dim(varSeqs_d2)
#3074271


################################################################################
#Determine Shared
################################################################################

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
shared_variant_matrix_withSexChroms <- shared_variant_matrix[!is.na(shared_variant_matrix[,1]),]

shared_variant_matrix_withSexChroms <- shared_variant_matrix_withSexChroms[!is.na(shared_variant_matrix_withSexChroms[,1]),]
dim(shared_variant_matrix_withSexChroms)
#809297

################################################################################
#Remove sex chromosomes.
################################################################################
toRemove <- grep("chrX", varSeqs_d1[,1])
if(length(toRemove)>0) {
  varSeqs_d1 <- varSeqs_d1[-toRemove,]
}

toRemove <- grep("chrX", varSeqs_d2[,1])
if(length(toRemove)>0) {
  varSeqs_d2 <- varSeqs_d2[-toRemove,]
}
toRemove <- grep("chrY", varSeqs_d1[,1])
if(length(toRemove)>0) {
  varSeqs_d1 <- varSeqs_d1[-toRemove,]
}

toRemove <- grep("chrY", varSeqs_d2[,1])
if(length(toRemove)>0) {
  varSeqs_d2 <- varSeqs_d2[-toRemove,]
}

dim(varSeqs_d1)
#3032858
dim(varSeqs_d2)
#3074271


################################################################################
#Determine Shared
################################################################################


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

dim(shared_variant_matrix)
#809297


################################################################################
#Restrict to heterozygous cases
################################################################################

varSeqs_d1_het <- varSeqs_d1[grep("HET", varSeqs_d1[,"HapM"]),]
varSeqs_d2_het <- varSeqs_d2[grep("HET", varSeqs_d2[,"HapM"]),]

dim(varSeqs_d1_het)
#1894277
dim(varSeqs_d2_het)
#1932258


################################################################################
#Restrict shared to heterozygous cases
################################################################################

varSeqs_d1_het_shared <- varSeqs_d1_het[varSeqs_d1_het[,"Region"]%in%varSeqs_d2_het[,"Region"],]
varSeqs_d2_het_shared <- varSeqs_d2_het[varSeqs_d2_het[,"Region"]%in%varSeqs_d1_het[,"Region"],]

varSeqs_d1_het_shared <- varSeqs_d1_het_shared[order(varSeqs_d1_het_shared[,"Region"]),]
varSeqs_d2_het_shared <- varSeqs_d2_het_shared[order(varSeqs_d2_het_shared[,"Region"]),]

shared_variant_matrix_hets <- matrix(nrow=nrow(varSeqs_d1_het_shared), ncol=6, data=NA)
for (i in 1:nrow(varSeqs_d1_het_shared)) {
  if(varSeqs_d1_het_shared[i,1]==varSeqs_d2_het_shared[i,1] && varSeqs_d1_het_shared[i,2]==varSeqs_d2_het_shared[i,2]) {
    shared_variant_matrix_hets[i,1] <- varSeqs_d1_het_shared[i,"Region"]
    shared_variant_matrix_hets[i,2] <- 1
    shared_variant_matrix_hets[i,3] <- varSeqs_d1_het_shared[i,"HapM"]
    shared_variant_matrix_hets[i,4] <- varSeqs_d1_het_shared[i,"HapP"]
    shared_variant_matrix_hets[i,5] <- varSeqs_d1_het_shared[i,"HapM"]
    shared_variant_matrix_hets[i,6] <- varSeqs_d1_het_shared[i,"HapP"]
  }
  if(varSeqs_d1_het_shared[i,1]==varSeqs_d2_het_shared[i,2] && varSeqs_d1_het_shared[i,2]==varSeqs_d2_het_shared[i,1]) {
    shared_variant_matrix_hets[i,1] <- varSeqs_d2_het_shared[i,"Region"]
    shared_variant_matrix_hets[i,2] <- 2
    shared_variant_matrix_hets[i,3] <- varSeqs_d2_het_shared[i,"HapM"]
    shared_variant_matrix_hets[i,4] <- varSeqs_d2_het_shared[i,"HapP"]
    shared_variant_matrix_hets[i,5] <- varSeqs_d2_het_shared[i,"HapM"]
    shared_variant_matrix_hets[i,6] <- varSeqs_d2_het_shared[i,"HapP"]
  }
}
shared_variant_matrix_hets <- shared_variant_matrix_hets[!is.na(shared_variant_matrix_hets[,1]),]

dim(shared_variant_matrix_hets)
#266448


################################################################################
#Read in pval tables summed across tissues. Remove inappropriate entries.
################################################################################

pval_table_1 <- read.table(pval_table_1_summed_saveFile, header=T, sep="\t", stringsAsFactors=F)

pval_table_2 <- read.table(pval_table_2_summed_saveFile, header=T, sep="\t", stringsAsFactors=F)

pval_table_1 <- pval_table_1[,!colnames(pval_table_1)=="INPUT"]
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)=="INPUT"]

toRemove <- grep("chrX", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
}

toRemove <- grep("chrX", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
}
toRemove <- grep("chrY", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
}

toRemove <- grep("chrY", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
}


################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.05 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#regardless of the factor.
#Do the same after removing cases of significant input.
################################################################################

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1))<=0.05]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2))<=0.05]

length(sig_d1)
#351094
length(sig_d2)
#205570

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#61429

sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#349406
length(sig_d2_noInput)
#205237

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#61276

################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.05 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#Excluding histones and POL2.
#Do the same after removing cases of significant input.
################################################################################

forbiddenExprs <- c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]))<=0.05]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs]))<=0.05]

length(sig_d1)
#136633
length(sig_d2)
#140709

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#28183



sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#135118
length(sig_d2_noInput)
#140399

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#28048


################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.001 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#regardless of the factor.
#Do the same after removing cases of significant input.
################################################################################

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1))<=0.001]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2))<=0.001]

length(sig_d1)
#86559
length(sig_d2)
#22074


dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#12121



sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#85245
length(sig_d2_noInput)
#21889

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#12020


################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.001 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#Excluding histones and POL2.
#Do the same after removing cases of significant input.
################################################################################

forbiddenExprs <- c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]))<=0.001]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs]))<=0.001]

length(sig_d1)
#8848
length(sig_d2)
#9733

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#1272

a <- sig_d1[sig_d1%in%sig_d2]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,])
#420




sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#7828
length(sig_d2_noInput)
#9570

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#1195


a <- sig_d1_noInput[sig_d1_noInput%in%sig_d2_noInput]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,])
#377


################################################################################
#Identify significant cases in histones and non-histones separately
#at P<=0.001 for each donor.
#Determine the number of variant significant for ONLY TFs in
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#Do the same after removing cases of significant input.
################################################################################


sig_d1_hists <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1[,colnames(pval_table_1)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")]))<=0.001]
sig_d2_hists <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2[,colnames(pval_table_2)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")]))<=0.001]


forbiddenExprs <- c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]))<=0.001]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs]))<=0.001]

sig_d1 <- sig_d1[!sig_d1%in%sig_d1_hists]
sig_d2 <- sig_d2[!sig_d2%in%sig_d2_hists]


length(sig_d1)
#4478
length(sig_d2)
#5809

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#904



sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#4448
length(sig_d2_noInput)
#5771

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#887




################################################################################
#Read in significance tables which are not summed across tissues, get
#number of sig variants (shared and unshared) in
################################################################################


pval_table_1 <- read.table(pval_table_1_all_saveFile, header=T, sep="\t", stringsAsFactors=F)

pval_table_2 <- read.table(pval_table_2_all_saveFile, header=T, sep="\t", stringsAsFactors=F)

pval_table_1 <- pval_table_1[,-grep("INPUT", colnames(pval_table_1))]
pval_table_2 <- pval_table_2[,-grep("INPUT", colnames(pval_table_2))]

toRemove <- grep("chrX", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
}

toRemove <- grep("chrX", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
}
toRemove <- grep("chrY", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
}

toRemove <- grep("chrY", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
}



################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.05 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#regardless of the factor.
#Do the same after removing cases of significant input.
################################################################################

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1))<=0.05]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2))<=0.05]

length(sig_d1)
#176397
length(sig_d2)
#166076

a <- shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),]
a <- a[!is.na(a[,1]),]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#35270


sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#174826
length(sig_d2_noInput)
#165757


dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#35126



################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.05 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#Excluding histones and POL2.
#Do the same after removing cases of significant input.
################################################################################

forbiddenExprs <- colnames(pval_table_1)[c(grep("H3K27AC", colnames(pval_table_1)), grep("H3K27ME3", colnames(pval_table_1)), grep("H3K4ME1", colnames(pval_table_1)), grep("H3K9AC", colnames(pval_table_1)), grep("POL2", colnames(pval_table_1)))]
forbiddenExprs_2 <- colnames(pval_table_2)[c(grep("H3K27AC", colnames(pval_table_2)), grep("H3K27ME3", colnames(pval_table_2)), grep("H3K4ME1", colnames(pval_table_2)), grep("H3K9AC", colnames(pval_table_2)), grep("POL2", colnames(pval_table_2)))]

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]))<=0.05]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs_2]))<=0.05]

length(sig_d1)
#140748
length(sig_d2)
#145261

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#28889

a <- sig_d1[sig_d1%in%sig_d2]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,])
#5954


sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#139234
length(sig_d2_noInput)
#144952


dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#28751

a <- sig_d1_noInput[sig_d1_noInput%in%sig_d2_noInput]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,])
#5862



################################################################################
#Determine the number and fraction of sig variants which are shared in each donor.
################################################################################

length(sig_d1[sig_d1%in%shared_variant_matrix_hets[,1]])
#17021
length(sig_d1[sig_d1%in%shared_variant_matrix_hets[,1]]) / length(sig_d1)
#0.1209324
#12.09%


length(sig_d2[sig_d2%in%shared_variant_matrix_hets[,1]])
#17822
length(sig_d2[sig_d2%in%shared_variant_matrix_hets[,1]]) / length(sig_d2)
#0.1226895
#12.26%



a <- sig_d1[sig_d1%in%sig_d2]
a <- a[a%in%shared_variant_matrix_hets[,1]]
length(a)
#5954



################################################################################
#Perform a bootstrapping test for this. Repeatedly subsample 28889 from the shared
#variants and determine the number of total significant variants found each time.
#10,000 samplings. Min pval of 0.0001.
################################################################################


sampleset_d1 <- varSeqs_d1_het[,1]
sampleset_d2 <- varSeqs_d2_het[,1]



length_vec <- c()
for (i in 1:10000) {
  if(i%%100==0) {print(i)}
  a <- sample(sampleset_d1, length(sig_d1), replace=F)
  b <- sample(sampleset_d2, length(sig_d2), replace=F)
  thisLength <- nrow(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(a, b),])
  length_vec <- c(length_vec, thisLength)
}
length(length_vec[length_vec<28889])/length(length_vec)


length_vec_alt <- c()
for (i in 1:10000) {
  if(i%%100==0) {print(i)}
  a <- sample(sampleset_d1, length(sig_d1), replace=F)
  b <- sample(sampleset_d2, length(sig_d2), replace=F)
  a_1 <- shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,1]
  b_1 <- shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%b,1]
  a_alt_number <- length(a[a%in%b])
  length_vec_alt <- c(length_vec_alt, a_alt_number)
}
length(length_vec_alt[length_vec_alt<alt_number])/length(length_vec_alt)










################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.001 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#regardless of the factor.
#Do the same after removing cases of significant input.
################################################################################

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1))<=0.001]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2))<=0.001]

length(sig_d1)
#5840
length(sig_d2)
#4565

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#628


sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#5196
length(sig_d2_noInput)
#4466


dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#575




################################################################################
#Determine for all experiments the number of variants which are significant
#at p<=0.001 and:
#1) Just donor1
#2) Just donor2
#3) Shared (significant in either case)
#Excluding histones and POL2.
#Do the same after removing cases of significant input.
################################################################################

forbiddenExprs <- colnames(pval_table_1)[c(grep("H3K27AC", colnames(pval_table_1)), grep("H3K27ME3", colnames(pval_table_1)), grep("H3K4ME1", colnames(pval_table_1)), grep("H3K9AC", colnames(pval_table_1)), grep("POL2", colnames(pval_table_1)))]
forbiddenExprs_2 <- colnames(pval_table_2)[c(grep("H3K27AC", colnames(pval_table_2)), grep("H3K27ME3", colnames(pval_table_2)), grep("H3K4ME1", colnames(pval_table_2)), grep("H3K9AC", colnames(pval_table_2)), grep("POL2", colnames(pval_table_2)))]

sig_d1 <- rownames(pval_table_1)[rowMins(as.matrix(pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]))<=0.001]
sig_d2 <- rownames(pval_table_2)[rowMins(as.matrix(pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs_2]))<=0.001]

length(sig_d1)
#4869
length(sig_d2)
#3970

dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1, sig_d2),])
#530


a <- sig_d1[sig_d1%in%sig_d2]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,])
#161



sig_d1_noInput <- sig_d1[!sig_d1%in%sigInput_1]
sig_d2_noInput <- sig_d2[!sig_d2%in%sigInput_1]

length(sig_d1_noInput)
#4328
length(sig_d2_noInput)
#3876


dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%c(sig_d1_noInput, sig_d2_noInput),])
#486


a <- sig_d1_noInput[sig_d1_noInput%in%sig_d2_noInput]
dim(shared_variant_matrix_hets[shared_variant_matrix_hets[,1]%in%a,])
#142



#
