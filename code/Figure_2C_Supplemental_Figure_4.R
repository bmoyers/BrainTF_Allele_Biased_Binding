#!/usr/bin/env R
#Figure_2C_Supplemental_Figure_4.R



################################################################################
#This script is used to produce Figures 2C and Supplemental 4.
#
#Run under R 4.1.0, but can be used under other versions so long as the packages
#     GenomicRanges, matrixStats, and ggplot are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# pval_table_1: A table compiling all p-values for each TF over all variants
#     in each dataset for donor 1, created by the script
#     CompilingHaplotypePvalsTable.R
# pval_table_2: As pval_table_1, but for donor 2
#
################################################################################


#module load cluster/R/4.1.0
#cd /cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/



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
pval_table_1 <- args[2]
pval_table_2 <- args[3]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg.txt"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg.txt"



################################################################################
#I want to show high correlation across regions for bias.
#Easiest way to do this seems to be p-value.
#So, for every TF and haplotype, determine if any tissue is significant at the
#1e-4 level.
#Do this for TFs done in only 4 tissues, as well as TFs done in all 9...
#Get the p-values for everything at that level.
#Get -log10 pval.
#make pairwise corr plot...
################################################################################






################################################################################
#Next, load in pval tables. Remove pol and H3K cases.
################################################################################

pval_table_1 <- read.table(pval_table_1, header=T, sep="\t", stringsAsFactors=F)
pval_table_1 <- pval_table_1[-grep("chrX", rownames(pval_table_1)),]

sigInput_d1 <- pval_table_1[which(rowMins(as.matrix(pval_table_1[,grep("INPUT", colnames(pval_table_1))]))<=0.05),]

allMat <- matrix(nrow=ncol(pval_table_1), ncol=3, data=unlist(strsplit(colnames(pval_table_1), split="_")), byrow=T)
allTFs <- unique(allMat[,2])
allTissues <- unique(allMat[,1])

forbiddenExprs <- c("CB_ZNF207_1224", "CB_ASCL1_1224", "CB_ATF7_1224", paste(allTissues, "OLIG2_1224", sep="_"), "CB_OLIG2_1230", "CB_ASCL1_1224")
if(length(forbiddenExprs[forbiddenExprs%in%colnames(pval_table_1)])>0) { pval_table_1 <- pval_table_1[,-which(colnames(pval_table_1)%in%forbiddenExprs)]}

forbiddenExprs_2 <- c(colnames(pval_table_1)[grep("H3K", colnames(pval_table_1))], colnames(pval_table_1)[grep("INPUT", colnames(pval_table_1))], colnames(pval_table_1)[grep("POL2", colnames(pval_table_1))])
pval_table_1 <- pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs_2]

set_4 <- c("CB", "DLPFC", "FP", "OL")


finalSet_d1 <- c()
finalSet_4_d1 <- c()

for (i in 1:length(allTFs)) {
  print(paste(i, length(allTFs)))
  thisSet <- pval_table_1[,grep(paste("_", allTFs[i], "_", sep=""), colnames(pval_table_1))]
  if(ncol(thisSet)>0) {
    #thisSet <- thisSet[rowMins(as.matrix(thisSet))<=0.05,]
    thisSet <- thisSet[rowMins(as.matrix(thisSet))<=0.001,]
    thisSet_Regions <- c()
    for (j in 1:ncol(thisSet)) { thisSet_Regions <- c(thisSet_Regions, strsplit(colnames(thisSet)[j], split="_")[[1]][1])}
    colnames(thisSet) <- thisSet_Regions
    for (j in 1:ncol(thisSet)) { thisSet[,j] <- -log10(thisSet[,j])}

    thisSet <- thisSet[,order(colnames(thisSet))]
    thisSet_4 <- thisSet[,colnames(thisSet)%in%set_4]

    if(ncol(thisSet)==9) { finalSet_d1 <- rbind(finalSet_d1, thisSet)}
    if(ncol(thisSet_4)==4) { finalSet_4_d1 <- rbind(finalSet_4_d1, thisSet_4)}
  }
}

finalSet_d1_noInput <- finalSet_d1[!rownames(finalSet_d1)%in%rownames(sigInput_d1),]
finalSet_4_d1_noInput <- finalSet_4_d1[!rownames(finalSet_4_d1)%in%rownames(sigInput_d1),]


cor(finalSet_d1)
cor(finalSet_4_d1)


cor(finalSet_d1_noInput)
cor(finalSet_4_d1_noInput)





pval_table_2 <- read.table(pval_table_2, header=T, sep="\t", stringsAsFactors=F)

sigInput_d2 <- pval_table_2[which(rowMins(as.matrix(pval_table_2[,grep("INPUT", colnames(pval_table_2))]))<=0.05),]

allMat <- matrix(nrow=ncol(pval_table_2), ncol=3, data=unlist(strsplit(colnames(pval_table_2), split="_")), byrow=T)
allTFs <- unique(allMat[,2])
allTissues <- unique(allMat[,1])

forbiddenExprs <- c("CB_ZNF207_1224", "CB_ASCL1_1224", "CB_ATF7_1224", paste(allTissues, "OLIG2_1224", sep="_"), "CB_OLIG2_1230", "CB_ASCL1_1224")
if(length(forbiddenExprs[forbiddenExprs%in%colnames(pval_table_2)])>0) { pval_table_2 <- pval_table_2[,-which(colnames(pval_table_2)%in%forbiddenExprs)]}

forbiddenExprs_2 <- c(colnames(pval_table_2)[grep("H3K", colnames(pval_table_2))], colnames(pval_table_2)[grep("INPUT", colnames(pval_table_2))], colnames(pval_table_2)[grep("POL2", colnames(pval_table_2))])
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs_2]

set_4 <- c("CB", "DLPFC", "FP", "OL")


finalSet_d2 <- c()
finalSet_4_d2 <- c()

for (i in 1:length(allTFs)) {
  print(paste(i, length(allTFs)))
  thisSet <- pval_table_2[,grep(paste("_", allTFs[i], "_", sep=""), colnames(pval_table_2))]
  if(ncol(thisSet)>0) {
    #thisSet <- thisSet[rowMins(as.matrix(thisSet))<=0.05,]
    thisSet <- thisSet[rowMins(as.matrix(thisSet))<=0.001,]
    thisSet_Regions <- c()
    for (j in 1:ncol(thisSet)) { thisSet_Regions <- c(thisSet_Regions, strsplit(colnames(thisSet)[j], split="_")[[1]][1])}
    colnames(thisSet) <- thisSet_Regions
    for (j in 1:ncol(thisSet)) { thisSet[,j] <- -log10(thisSet[,j])}

    thisSet <- thisSet[,order(colnames(thisSet))]
    thisSet_4 <- thisSet[,colnames(thisSet)%in%set_4]

    if(ncol(thisSet)==9) { finalSet_d2 <- rbind(finalSet_d2, thisSet)}
    if(ncol(thisSet_4)==4) { finalSet_4_d2 <- rbind(finalSet_4_d2, thisSet_4)}
  }
}


finalSet_d2_noInput <- finalSet_d2[!rownames(finalSet_d2)%in%rownames(sigInput_d2),]
finalSet_4_d2_noInput <- finalSet_4_d2[!rownames(finalSet_4_d2)%in%rownames(sigInput_d2),]


cor(finalSet_d2)
cor(finalSet_4_d2)

cor(finalSet_d2_noInput)
cor(finalSet_4_d2_noInput)




finalSet <- rbind(finalSet_d1, finalSet_d2)
finalSet_4 <- rbind(finalSet_4_d1, finalSet_4_d2)


finalSet_noInput <- rbind(finalSet_d1_noInput, finalSet_d2_noInput)
finalSet_4_noInput <- rbind(finalSet_4_d1_noInput, finalSet_4_d2_noInput)


cor(finalSet)
cor(finalSet_4)

cor(finalSet_noInput)
cor(finalSet_4_noInput)




panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- abs(cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}


saveFile <- paste(outDir, "Figure_2C.pdf", sep="")
pdf(saveFile)
pairs(finalSet_noInput, upper.panel = panel.cor, lower.panel = panel.smooth)
dev.off()


saveFile <- paste(outDir, "Supplemental_Figure_3.pdf", sep="")
pdf(saveFile)
pairs(finalSet_4_noInput, upper.panel = panel.cor, lower.panel = panel.smooth)
dev.off()








#
