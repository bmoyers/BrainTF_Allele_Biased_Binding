#!/usr/bin/env R
#Supplemental_5_6.R


################################################################################
#This script is used to produce Supplemental Figures 4 and 5.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges and matrixStats packages are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# pvals_1: A table compiling p-value calculations when reads are summed across
#     experiments for donor1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# pvals_2: As pval_table_1, but for donor 2.
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(matrixStats)
library(ggplot2)

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
pvals_1 <- args[2]
pvals_2 <- args[3]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#pvals_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pvals_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"

################################################################################
#Load in the p-value tables.
################################################################################

pvals_1 <- read.table(pvals_1, header=T, sep="\t", stringsAsFactors=F)

pvals_2 <- read.table(pvals_2, header=T, sep="\t", stringsAsFactors=F)


################################################################################
#Remove all cases of the x-chromosome.
################################################################################

if(length(grep("chrX", rownames(pvals_1)))>0) {
  pvals_1 <- pvals_1[-grep("chrX", rownames(pvals_1)),]
}

if(length(grep("chrX", rownames(pvals_2)))>0) {
  pvals_2 <- pvals_2[-grep("chrX", rownames(pvals_2)),]
}


################################################################################
#Identify cases where Input was significant, and remove them.
################################################################################

sigInput_1 <- which(pvals_1$INPUT<=0.05)
if(length(sigInput_1)>0) { pvals_1 <- pvals_1[-sigInput_1,] }


sigInput_2 <- which(pvals_2$INPUT<=0.05)
if(length(sigInput_2)>0) { pvals_2 <- pvals_2[-sigInput_2,] }



################################################################################
#Exclude input.
################################################################################

pvals_1 <- pvals_1[,!colnames(pvals_1)%in%c("INPUT")]
pvals_2 <- pvals_2[,!colnames(pvals_2)%in%c("INPUT")]


################################################################################
#Identify cases where at least one factor is significant. Restrict to those.
################################################################################

pvals_1 <- pvals_1[rowMins(as.matrix(pvals_1))<=0.001,]
pvals_2 <- pvals_2[rowMins(as.matrix(pvals_2))<=0.001,]


################################################################################
#Identify cases where TFs are significant.
################################################################################

sig_tf_1 <- which(rowMins(as.matrix(pvals_1[,!colnames(pvals_1)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")]))<=0.001)
sig_tf_2 <- which(rowMins(as.matrix(pvals_2[,!colnames(pvals_2)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")]))<=0.001)


################################################################################
#Identify cases where histones are significant.
################################################################################

sig_hist_1 <- which(rowMins(as.matrix(pvals_1[,colnames(pvals_1)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")]))<=0.001)
sig_hist_2 <- which(rowMins(as.matrix(pvals_2[,colnames(pvals_2)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "POL2")]))<=0.001)


################################################################################
#For cases where histones ARE NOT significant, identify:
#1) the number of TFs which are significant
#2) all p-values (-log10?) of significant experiments.
################################################################################

sig_tf_noHist_1 <- sig_tf_1[!sig_tf_1%in%sig_hist_1]
sig_tf_noHist_2 <- sig_tf_2[!sig_tf_2%in%sig_hist_2]

numSig_TFonly_1 <- c()
pval_TFonly_1 <- c()
for (i in 1:length(sig_tf_noHist_1)) {
  thisSet <- as.numeric(pvals_1[sig_tf_noHist_1[i],!colnames(pvals_1)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")])
  numSig_TFonly_1 <- c(numSig_TFonly_1, length(thisSet[thisSet<=0.001]))
  pval_TFonly_1 <- c(pval_TFonly_1, -log10(thisSet[thisSet<=0.001]+0.000000000000000000001))
}

numSig_TFonly_2 <- c()
pval_TFonly_2 <- c()
for (i in 1:length(sig_tf_noHist_2)) {
  thisSet <- as.numeric(pvals_2[sig_tf_noHist_2[i],!colnames(pvals_2)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")])
  numSig_TFonly_2 <- c(numSig_TFonly_2, length(thisSet[thisSet<=0.001]))
  pval_TFonly_2 <- c(pval_TFonly_2, -log10(thisSet[thisSet<=0.001]+0.000000000000000000001))
}


################################################################################
#For cases where histones ARE significant, identify:
#1) the number of TFs which are significant
#2) all p-values (-log10?) of significant TF experiments.
################################################################################

sig_tf_Hist_1 <- sig_tf_1[sig_tf_1%in%sig_hist_1]
sig_tf_Hist_2 <- sig_tf_2[sig_tf_2%in%sig_hist_2]

numSig_TFHist_1 <- c()
pval_TFHist_1 <- c()
for (i in 1:length(sig_tf_Hist_1)) {
  thisSet <- as.numeric(pvals_1[sig_tf_Hist_1[i],!colnames(pvals_1)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")])
  numSig_TFHist_1 <- c(numSig_TFHist_1, length(thisSet[thisSet<=0.001]))
  pval_TFHist_1 <- c(pval_TFHist_1, -log10(thisSet[thisSet<=0.001]+0.000000000000000000001))
}

numSig_TFHist_2 <- c()
pval_TFHist_2 <- c()
for (i in 1:length(sig_tf_Hist_2)) {
  thisSet <- as.numeric(pvals_2[sig_tf_Hist_2[i],!colnames(pvals_2)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")])
  numSig_TFHist_2 <- c(numSig_TFHist_2, length(thisSet[thisSet<=0.001]))
  pval_TFHist_2 <- c(pval_TFHist_2, -log10(thisSet[thisSet<=0.001]+0.000000000000000000001))
}


################################################################################
#Make relevant density plots.
################################################################################

numSig_TFs <- c(numSig_TFonly_1, numSig_TFonly_2)
numSig_Hists <- c(numSig_TFHist_1, numSig_TFHist_2)

numSig_TFs <- cbind(numSig_TFs, rep("TF Only", length(numSig_TFs)))
colnames(numSig_TFs) <- c("number_sig", "type")

numSig_Hists <- cbind(numSig_Hists, rep("TF with Hist", length(numSig_Hists)))
colnames(numSig_Hists) <- c("number_sig", "type")

graphDF <- rbind(numSig_TFs, numSig_Hists)
graphDF <- as.data.frame(graphDF)
graphDF[,1] <- as.numeric(graphDF[,1])

saveFile <- paste(outDir, "Supplemental_Figure_4.pdf", sep="")
tfNumDensity <- ggplot(graphDF, aes(number_sig, fill=type)) +
  geom_bar(position="dodge") + xlim(0,20) +
  theme_classic(base_size=20)
ggsave(saveFile)




pvalSig_TFs <- c(pval_TFonly_1, pval_TFonly_2)
pvalSig_Hists <- c(pval_TFHist_1, pval_TFHist_2)

pvalSig_TFs <- cbind(pvalSig_TFs, rep("TF Only", length(pvalSig_TFs)))
colnames(pvalSig_TFs) <- c("neg_log10_pval", "type")

pvalSig_Hists <- cbind(pvalSig_Hists, rep("TF with Hist", length(pvalSig_Hists)))
colnames(pvalSig_Hists) <- c("neg_log10_pval", "type")

graphDF <- rbind(pvalSig_TFs, pvalSig_Hists)
graphDF <- as.data.frame(graphDF)
graphDF[,1] <- as.numeric(graphDF[,1])

saveFile <- paste(outDir, "Supplemental_Figure_5.pdf", sep="")
tfNumDensity <- ggplot(graphDF, aes(neg_log10_pval, fill=type)) +
  geom_density(alpha=0.5) +
  theme_classic(base_size=20)
ggsave(saveFile)


################################################################################
#Related, but separate analysis-- for cases which are significant for
#at least one TF, determine the percentage that is significant for
#multiple TFs.
#Determine also the percentage which are significant for a histone or POL.
################################################################################

count_multiTF_d1 <- 0
count_histEffect_d1 <- 0

for (i in 1:length(sig_tf_1)) {
  if(sig_tf_1[i]%in%sig_hist_1) { count_histEffect_d1 <- count_histEffect_d1+1 }
  thisSet_tfs_only <- pvals_1[sig_tf_1[i],!colnames(pvals_1)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")]
  thisSet_tfs_only <- as.numeric(thisSet_tfs_only)
  thisCount <- length(thisSet_tfs_only[thisSet_tfs_only<=0.001])
  if(thisCount>=2) { count_multiTF_d1 <- count_multiTF_d1+1}
}



count_multiTF_d2 <- 0
count_histEffect_d2 <- 0

for (i in 1:length(sig_tf_2)) {
  if(sig_tf_2[i]%in%sig_tf_2) { count_histEffect_d2 <- count_histEffect_d2+1 }
  thisSet_tfs_only <- pvals_2[sig_tf_2[i],!colnames(pvals_2)%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC")]
  thisSet_tfs_only <- as.numeric(thisSet_tfs_only)
  thisCount <- length(thisSet_tfs_only[thisSet_tfs_only<=0.001])
  if(thisCount>=2) { count_multiTF_d2 <- count_multiTF_d2+1}
}


multi_tf_perc <- (count_multiTF_d1 + count_multiTF_d2) / (length(sig_tf_1) + length(sig_tf_2))
tf_hist_perc <- (count_histEffect_d1 + count_histEffect_d2) / (length(sig_tf_1) + length(sig_tf_2))




#
