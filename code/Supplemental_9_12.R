#!/usr/bin/env R
#Supplemental_9_12.R


################################################################################
#This script is used to produce Supplemental Figures 9 and 12.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges, ggplot2, and matrixStats packages are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# pval_table_1: A table compiling p-value calculations across all experiments for
#     donor 1 when summed across tissues.
# pval_table_2: A table compiling p-value calculations across all experiments for
#     donor 2 when summed across tissues.
# theMotifs: The JASPAR motif database.
#
################################################################################


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(matrixStats)
library(ggplot2)
library(GenomicRanges)
library(memes)
library(universalmotif)

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
pval_table_1 <- args[2]
pval_table_2 <- args[3]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#theMotifs <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/jaspar_2023July03/JASPAR2022_CORE_non-redundant_pfms_meme.txt"




################################################################################
#Load in p-value tables; Remove cases of significant input.  Remove
#cases of significant input.  Remove Histones and Pol.
################################################################################


pval_table_1 <- read.delim(pval_table_1, header=T, sep="\t", stringsAsFactors=F)
pval_table_2 <- read.delim(pval_table_2, header=T, sep="\t", stringsAsFactors=F)

sigInput_1 <- rownames(pval_table_1[which(pval_table_1$INPUT<=0.05),])
sigInput_2 <- rownames(pval_table_2[which(pval_table_2$INPUT<=0.05),])

pval_table_1 <- pval_table_1[!rownames(pval_table_1)%in%sigInput_1,]
pval_table_2 <- pval_table_2[!rownames(pval_table_2)%in%sigInput_2,]

forbiddenExprs_2 <- c(colnames(pval_table_1)[grep("H[3|4]K", colnames(pval_table_1))], colnames(pval_table_1)[grep("INPUT", colnames(pval_table_1))], colnames(pval_table_1)[grep("POL2", colnames(pval_table_1))])
pval_table_1 <- pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs_2]

forbiddenExprs_2 <- c(colnames(pval_table_2)[grep("H[3|4]K", colnames(pval_table_2))], colnames(pval_table_2)[grep("INPUT", colnames(pval_table_2))], colnames(pval_table_2)[grep("POL2", colnames(pval_table_2))])
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs_2]


pval_table_1 <- pval_table_1[rowMins(as.matrix(pval_table_1))<=0.001,]
pval_table_2 <- pval_table_2[rowMins(as.matrix(pval_table_2))<=0.001,]




################################################################################
#Count the total number of significant regions for each donor.
#For each TF, count the number of cases where that TF is significant.
################################################################################

uniqueTFs <- unique(colnames(pval_table_2))

fractionTable <- c()
for (i in 1:length(uniqueTFs)) {
  if(uniqueTFs[i]%in%colnames(pval_table_1)) {
    firstNum <- nrow(pval_table_1[pval_table_1[,uniqueTFs[i]]<=0.001,])
    secondNum <- nrow(pval_table_2[pval_table_2[,uniqueTFs[i]]<=0.001,])
    totalFraction <- (firstNum + secondNum) / (nrow(pval_table_1) + nrow(pval_table_2))
    fractionTable <- rbind(fractionTable, c(uniqueTFs[i], firstNum, nrow(pval_table_1), secondNum, nrow(pval_table_2), totalFraction))
  }
  if(!uniqueTFs[i]%in%colnames(pval_table_1)) {
    firstNum <- NA
    secondNum <- nrow(pval_table_2[pval_table_2[,uniqueTFs[i]]<=0.001,])
    totalFraction <- (secondNum) / (nrow(pval_table_2))
    fractionTable <- rbind(fractionTable, c(uniqueTFs[i], firstNum, NA, secondNum, nrow(pval_table_2), totalFraction))
  }

}

fractionTable <- as.data.frame(fractionTable)
colnames(fractionTable) <- c("TF", "TF_Sig_d1", "Total_Sig_d1", "TF_Sig_d2", "Total_sig_d2", "Fraction")

fractionTable$Fraction <- as.numeric(as.character(fractionTable$Fraction))
fractionTable <- fractionTable[order(fractionTable$Fraction, decreasing=T),]
fractionTable$TF <- factor(fractionTable$TF, levels=unique(fractionTable$TF))



################################################################################
#Read in the motifs.  Add information about whether or not a motif for the TF is
#present in the database to color the final plot. After examining it, this
#appears to not be productive, so removing the coloring.
################################################################################

jaspar_motifs <- read_meme(theMotifs)
alt_jaspar_motifs <- filter_motifs(jaspar_motifs, altname=uniqueTFs)

in_jaspar <- c()
for (i in 1:length(alt_jaspar_motifs)) {
  in_jaspar <- c(in_jaspar, alt_jaspar_motifs[[i]]["altname"])
}
in_jaspar <- unique(in_jaspar)

fractionTable$Motif <- rep("No", nrow(fractionTable))
fractionTable[fractionTable$TF%in%in_jaspar,"Motif"] <- "Yes"

################################################################################
#Make the plot
################################################################################


saveFile <- paste(outDir, "Supplemental_9.pdf", sep="")
ggplot(fractionTable, aes(y=Fraction, x=TF)) +
  geom_bar(position="stack", stat="identity") +
  ggtitle("") + xlab("TF") + ylab("Fraction of Sig Regions Found Sig") +
  theme_classic(base_size=22) + theme(axis.text.x = element_text(angle = 90, size=10))
ggsave(saveFile, width=14)



################################################################################
################################################################################
#Read in supplemental Table 5.  For this set of variants, determine
#the same information as the plot above.
################################################################################
################################################################################

suppTable5 <- paste(outDir, "Supplemental_Table_5.txt", sep="")
suppTable5 <- read.table(suppTable5, header=T, sep="\t", stringsAsFactors=F)




firstSet <- suppTable5[suppTable5$d1_tf_significance=="Significant",]
firstSet_gr <- makeGRangesFromDataFrame(firstSet)


pval_table_1_regions <- as.data.frame(matrix(nrow=nrow(pval_table_1), ncol=3, data=unlist(strsplit(rownames(pval_table_1), split=":|-")), byrow=T))
colnames(pval_table_1_regions) <- c("chr", "start", "end")
pval_table_1_regions_gr <- makeGRangesFromDataFrame(pval_table_1_regions)

theIntersect <- as.data.frame(findOverlaps(pval_table_1_regions_gr, firstSet_gr))

pval_table_1_gtex <- pval_table_1[unique(theIntersect[,1]),]



secondSet <- suppTable5[suppTable5$d2_state=="Heterozygous",]
secondSet <- secondSet[secondSet$d2_tf_significance=="Significant",]
secondSet_gr <- makeGRangesFromDataFrame(secondSet)

pval_table_2_regions <- as.data.frame(matrix(nrow=nrow(pval_table_2), ncol=3, data=unlist(strsplit(rownames(pval_table_2), split=":|-")), byrow=T))
colnames(pval_table_2_regions) <- c("chr", "start", "end")
pval_table_2_regions_gr <- makeGRangesFromDataFrame(pval_table_2_regions)

theIntersect <- as.data.frame(findOverlaps(pval_table_2_regions_gr, secondSet_gr))

pval_table_2_gtex <- pval_table_2[unique(theIntersect[,1]),]





################################################################################
#Count the total number of significant regions for each donor.
#For each TF, count the number of cases where that TF is significant.
################################################################################


uniqueTFs <- unique(colnames(pval_table_2_gtex))

fractionTable_gtex <- c()
for (i in 1:length(uniqueTFs)) {
  if(uniqueTFs[i]%in%colnames(pval_table_1_gtex)) {
    firstNum <- nrow(pval_table_1_gtex[pval_table_1_gtex[,uniqueTFs[i]]<=0.001,])
    secondNum <- nrow(pval_table_2_gtex[pval_table_2_gtex[,uniqueTFs[i]]<=0.001,])
    totalFraction <- (firstNum + secondNum) / (nrow(pval_table_1_gtex) + nrow(pval_table_2_gtex))
    fractionTable_gtex <- rbind(fractionTable_gtex, c(uniqueTFs[i], firstNum, nrow(pval_table_1_gtex), secondNum, nrow(pval_table_2_gtex), totalFraction))
  }
  if(!uniqueTFs[i]%in%colnames(pval_table_1_gtex)) {
    firstNum <- NA
    secondNum <- nrow(pval_table_2_gtex[pval_table_2_gtex[,uniqueTFs[i]]<=0.001,])
    totalFraction <- (secondNum) / (nrow(pval_table_2_gtex))
    fractionTable_gtex <- rbind(fractionTable_gtex, c(uniqueTFs[i], firstNum, NA, secondNum, nrow(pval_table_2_gtex), totalFraction))
  }

}

fractionTable_gtex <- as.data.frame(fractionTable_gtex)
colnames(fractionTable_gtex) <- c("TF", "TF_Sig_d1", "Total_Sig_d1", "TF_Sig_d2", "Total_sig_d2", "Fraction")

fractionTable_gtex$Fraction <- as.numeric(as.character(fractionTable_gtex$Fraction))
fractionTable_gtex <- fractionTable_gtex[order(fractionTable_gtex$Fraction, decreasing=T),]
fractionTable_gtex$TF <- factor(fractionTable_gtex$TF, levels=unique(fractionTable_gtex$TF))



################################################################################
#Read in the motifs.  Add information about whether or not a motif for the TF is
#present in the database to color the final plot. After examining it, this
#appears to not be productive, so removing the coloring.
################################################################################

jaspar_motifs <- read_meme(theMotifs)
alt_jaspar_motifs <- filter_motifs(jaspar_motifs, altname=uniqueTFs)

in_jaspar <- c()
for (i in 1:length(alt_jaspar_motifs)) {
  in_jaspar <- c(in_jaspar, alt_jaspar_motifs[[i]]["altname"])
}
in_jaspar <- unique(in_jaspar)

fractionTable_gtex$Motif <- rep("No", nrow(fractionTable_gtex))
fractionTable_gtex[fractionTable_gtex$TF%in%in_jaspar,"Motif"] <- "Yes"

################################################################################
#Make the plot
################################################################################


saveFile <- paste(outDir, "Supplemental_12.pdf", sep="")
ggplot(fractionTable_gtex, aes(y=Fraction, x=TF)) +
  geom_bar(position="stack", stat="identity") +
  ggtitle("") + xlab("TF") + ylab("Fraction of Sig Regions Found Sig") +
  theme_classic(base_size=22) + theme(axis.text.x = element_text(angle = 90, size=10))
ggsave(saveFile, width=14)





#
