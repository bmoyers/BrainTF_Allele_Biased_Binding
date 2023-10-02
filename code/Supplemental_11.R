#Supplemental_11.R


################################################################################
#This script is used to make supplemental figure 11.
#
#Run under R 4.1.0, but can be used under other versions.
#
# outDir: Directory to which figures should be saved.
# pval_table_summed_1: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# pval_table_summed_2: As pval_table_summed_1 but for donor 2.
# cadd_annotations_1_saveFile: a vep-annotated vcf file for donor 1, produced by vep using the
#     vep108.ini file.
# cadd_annotations_2_saveFile: a vep-annotated vcf file for donor 2, produced by vep using the
#     vep108.ini file.
# depths_m_d1: A table compiling all TF read  depths for donor 1, haplotype 1,
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depths_p_d1: As depths_m_d1, but for haplotype 2.
# depths_m_d2: A table compiling all TF read  depths for donor 2, haplotype 1,
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depths_p_d2: As depths_m_d2, but for haplotype 2.
#
################################################################################


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################


library(GenomicRanges)
library(ggplot2)
library(viridis)
library(matrixStats)



################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

#Extracts relevant information from cadd table and formats it for easier
#calculation of enrichments downstream.
format_cadd_info <- function(cadd_annotations) {
  theColnames <- as.character(cadd_annotations[2,])
  cadd_annotations <- cadd_annotations[3:nrow(cadd_annotations),]
  colnames(cadd_annotations) <- theColnames
  cadd_annotations_reduced <- cadd_annotations[,c("#Chrom", "Pos", "Ref", "Alt", "GerpS", "RawScore", "PHRED")]
  cadd_annotations_reduced <- unique(cadd_annotations_reduced)

  cadd_annotations_reduced <- cbind(paste("chr", cadd_annotations_reduced[,1], sep=""), cadd_annotations_reduced[,2], cadd_annotations_reduced[,2], cadd_annotations_reduced[,3:7])
  colnames(cadd_annotations_reduced) <- c("chr", "start", "end", "ref", "alt", "GerpS", "RawScore", "PHRED")
  cadd_annotations_reduced <- as.data.frame(cadd_annotations_reduced)
  cadd_annotations_reduced[,1] <- as.character(cadd_annotations_reduced[,1])
  cadd_annotations_reduced[,4] <- as.character(cadd_annotations_reduced[,4])
  cadd_annotations_reduced[,5] <- as.character(cadd_annotations_reduced[,5])
  for (i in c(2,3,6,7,8)) { cadd_annotations_reduced[,i] <- as.numeric(as.character(cadd_annotations_reduced[,i]))}

  return(cadd_annotations_reduced)
}



#This function calculates enrichment for each of the relevant categories.
getEnrichmentTable_combined <- function(pvalTable_1, depthTable_1, enrichmentOverlap_1, enrichments_1, pvalTable_2, depthTable_2, enrichmentOverlap_2, enrichments_2, suffDepth=11) {
  #pvalTable_1 <- pval_table_summed_1
  #depthTable_1 <- depth_table_summed_1
  #enrichmentOverlap_1 <- conservation_overlaps_d1
  #enrichments_1 <- cadd_annotations_1
  #pvalTable_2 <- pval_table_summed_2
  #depthTable_2 <- depth_table_summed_2
  #enrichmentOverlap_2 <- conservation_overlaps_d2
  #enrichments_2 <- cadd_annotations_2
  #suffDepth <- 11

  retTable <- c()
  for (i in 1:ncol(pvalTable_1)) {

    considered_1 <- which(depthTable_1[,i]>=suffDepth)
    significants_1 <- which(pvalTable_1[,i]<=0.05)
    significants_1 <- significants_1[significants_1%in%considered_1]

    enrichments_considered_1 <- enrichments_1[enrichmentOverlap_1[considered_1,2],]
    enrichments_significant_1 <- enrichments_1[enrichmentOverlap_1[significants_1,2],]

    considered_2 <- which(depthTable_2[,i]>=suffDepth)
    significants_2 <- which(pvalTable_2[,i]<=0.05)
    significants_2 <- significants_2[significants_2%in%considered_2]

    enrichments_considered_2 <- enrichments_2[enrichmentOverlap_2[considered_2,2],]
    enrichments_significant_2 <- enrichments_2[enrichmentOverlap_2[significants_2,2],]

    considered <- c(considered_1, considered_2)
    significants <- c(significants_1, significants_2)
    enrichments_considered <- rbind(enrichments_considered_1, enrichments_considered_2)
    enrichments_significant <- rbind(enrichments_significant_1, enrichments_significant_2)

    considered_posGERP_num <- nrow(enrichments_considered[enrichments_considered[,"GerpRS"]>4,])
    significant_posGERP_num <- nrow(enrichments_significant[enrichments_significant[,"GerpRS"]>4,])
    posGERP_enrichment <- (significant_posGERP_num/nrow(enrichments_significant))/(considered_posGERP_num/nrow(enrichments_considered))
    #gerp_pval_1 <- chisq.test(x=c(significant_posGERP_num, considered_posGERP_num-significant_posGERP_num), y=c(nrow(enrichments_significant), nrow(enrichments_considered)-nrow(enrichments_significant)))
    gerp_pval_1 <- chisq.test(x=c(significant_posGERP_num, nrow(enrichments_significant)-significant_posGERP_num), y=c(considered_posGERP_num, nrow(enrichments_considered)-considered_posGERP_num))
    gerp_pval_1 <- gerp_pval_1$p.value


    considered_posCADD_num <- nrow(enrichments_considered[enrichments_considered[,"CADD_RAW"]>1,])
    significant_posCADD_num <- nrow(enrichments_significant[enrichments_significant[,"CADD_RAW"]>1,])
    posCADD_enrichment <- (significant_posCADD_num/nrow(enrichments_significant))/(considered_posCADD_num/nrow(enrichments_considered))
    #cadd_pval_1 <- chisq.test(x=c(significant_posCADD_num, considered_posCADD_num-significant_posCADD_num), y=c(nrow(enrichments_significant), nrow(enrichments_considered)-nrow(enrichments_significant)))
    cadd_pval_1 <- chisq.test(x=c(significant_posCADD_num, nrow(enrichments_significant)-significant_posCADD_num), y=c(considered_posCADD_num, nrow(enrichments_considered)-considered_posCADD_num))
    cadd_pval_1 <- cadd_pval_1$p.value


    thisVec <- c(colnames(pvalTable_1)[i], length(considered), length(significants),
      posGERP_enrichment, posCADD_enrichment, gerp_pval_1, cadd_pval_1)

    retTable <- rbind(retTable, thisVec)
  }

  colnames(retTable) <- c("TF", "possibleSides", "significantSites", "posGerpEnrichment", "posCaddEnrichment", "gp1", "cp1")
  rownames(retTable) <- c()
  retTable <- as.data.frame(retTable)
  retTable[,1] <- as.character(retTable[,1])
  for (i in 2:ncol(retTable)) { retTable[,i] <- as.numeric(as.character(retTable[,i]))}
  return(retTable)
}



#Plots enrichment table.
plot_enrichments <- function(enrichmentTable, saveBase) {
  enrichmentTable <- as.data.frame(enrichmentTable)
  enrichmentTable$TF <- paste(enrichmentTable$TF, " (", enrichmentTable$significantSites, " sig. of ", enrichmentTable$possibleSides, ")", sep="")

  #colnames(enrichmentTable)[4:5] <- c("GERP>4", "CADD>1")
  #
  for (i in 4:ncol(enrichmentTable)) {
    saveFile <- paste(saveBase, "_", colnames(enrichmentTable)[i], ".pdf", sep="")

    this_enrichmentTable <- enrichmentTable
    this_enrichmentTable[is.na(this_enrichmentTable[,i]),i] <- 0
    bestOrdering <- this_enrichmentTable
    bestOrdering <- bestOrdering[order(bestOrdering[,i]),]
    bestOrdering <- as.character(bestOrdering$TF)
    this_enrichmentTable$TF <- factor(this_enrichmentTable$TF, levels=bestOrdering)
    colnames(this_enrichmentTable)[i] <- "EnrFactor"
    this_enrichmentTable[i] <- log(this_enrichmentTable[i])

    if(i==4) {ylab <- "GERP>4 log(Enrichment)"}
    if(i==5) {ylab <- "CADD>1 log(Enrichment)"}

    ggplot(this_enrichmentTable, aes(y=EnrFactor, x=TF)) +
      geom_bar(position="stack", stat="identity") +
      scale_fill_viridis(discrete = T) +
      ggtitle("") + xlab("") + ylab(ylab) +
      theme_classic(base_size=22) + theme(axis.text.x = element_text(angle = 90, size=10))
      #theme(axis.text=element_text(size=20), axis.title=element_text(size=24,face="bold"), axis.ticks.x=element_blank())
    ggsave(saveFile, width=14)
  }
}



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################


##################################################################
#Identify the minimum necessary read count for a p-value of 0.001
#or lower.
##################################################################

###Indentify minimum depth for pval of 0.001 or lower.

foundMin <- FALSE
minDepth <- 1
while(!foundMin) {
  thisTest <- binom.test(0,minDepth, alternative="two.sided")
  if(thisTest$p.value<=0.001) { foundMin <- TRUE}
  if(thisTest$p.value>0.001) { minDepth <- minDepth+1}
}
#minDepth is 11.


#


##################################################################
#Identify the outDir
##################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
pval_table_summed_1 <- args[2]
pval_table_summed_2 <- args[3]
cadd_annotations_1 <- args[4]
cadd_annotations_2 <- args[5]
depths_m_d1 <- args[6]
depths_p_d1 <- args[7]
depths_m_d2 <- args[8]
depths_p_d2 <- args [9]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#pval_table_summed_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pval_table_summed_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#cadd_annotations_1_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0001_annotated.tsv.gz"
#cadd_annotations_2_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0002_annotated.tsv.gz"
#depths_m_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depths_p_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#depths_m_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depths_p_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"



##################################################################
#Read in and prepare data.
##################################################################


pval_table_summed_1 <- read.table(pval_table_summed_1, header=T, sep="\t", stringsAsFactors=F)
pval_table_summed_2 <- read.table(pval_table_summed_2, header=T, sep="\t", stringsAsFactors=F)

depths_m_d1 <- read.table(depths_m_d1)
depths_p_d1 <- read.table(depths_p_d1)
depth_table_summed_1 <- depths_m_d1 + depths_p_d1
rm(depths_m_d1, depths_p_d1)

depths_m_d2 <- read.table(depths_m_d2)
depths_p_d2 <- read.table(depths_p_d2)
depth_table_summed_2 <- depths_m_d2 + depths_p_d2
rm(depths_m_d2, depths_p_d2)

sigInput_1 <- rownames(pval_table_summed_1[which(pval_table_summed_1$INPUT<=0.05),])
sigInput_2 <- rownames(pval_table_summed_2[which(pval_table_summed_2$INPUT<=0.05),])

pval_table_summed_1 <- pval_table_summed_1[!rownames(pval_table_summed_1)%in%sigInput_1,]
pval_table_summed_2 <- pval_table_summed_2[!rownames(pval_table_summed_2)%in%sigInput_2,]

depth_table_summed_1 <- depth_table_summed_1[!rownames(depth_table_summed_1)%in%sigInput_1,]
depth_table_summed_2 <- depth_table_summed_2[!rownames(depth_table_summed_2)%in%sigInput_2,]


##################################################################
#CADD and other annotations.
##################################################################

cadd_annotations_1 <- readLines(cadd_annotations_1_saveFile)
cadd_annotations_1 <- cadd_annotations_1[-grep("^##", cadd_annotations_1)]
cadd_annotations_1_header <- cadd_annotations_1[grep("^#", cadd_annotations_1)]
cadd_annotations_1_header <- gsub("^#", "", cadd_annotations_1_header)
cadd_annotations_1_header <- strsplit(cadd_annotations_1_header, split="\t")[[1]]

cadd_annotations_1 <- read.delim(cadd_annotations_1_saveFile, header=F, sep="\t", stringsAsFactors=F, comment.char="#")
colnames(cadd_annotations_1) <- cadd_annotations_1_header

#
cadd_annotations_1 <- cadd_annotations_1[,c("Location", "GerpRS", "CADD_RAW", "CADD_PHRED")]

cadd_annotations_1_locations <- as.data.frame(matrix(nrow=nrow(cadd_annotations_1), ncol=2, data=unlist(strsplit(as.character(cadd_annotations_1$Location), split=":", fixed=T)), byrow=T))
colnames(cadd_annotations_1_locations) <- c("chr", "start")
withDash <- grep("-", cadd_annotations_1_locations$start)
cadd_annotations_1_locations[-withDash, "start"] <- paste(cadd_annotations_1_locations[-withDash, "start"], cadd_annotations_1_locations[-withDash, "start"], sep="-")
corrected <- as.data.frame(matrix(nrow=nrow(cadd_annotations_1_locations), ncol=2, data=unlist(strsplit(cadd_annotations_1_locations[,2], split="-")), byrow=T))
cadd_annotations_1_locations <- cbind(cadd_annotations_1_locations$chr, corrected)
colnames(cadd_annotations_1_locations) <- c("chr", "start", "end")

cadd_annotations_1 <- cbind(cadd_annotations_1_locations, cadd_annotations_1[,2:4])
cadd_annotations_1_gr <- makeGRangesFromDataFrame(cadd_annotations_1)



cadd_annotations_2 <- readLines(cadd_annotations_2_saveFile)
cadd_annotations_2 <- cadd_annotations_2[-grep("^##", cadd_annotations_2)]
cadd_annotations_2_header <- cadd_annotations_2[grep("^#", cadd_annotations_2)]
cadd_annotations_2_header <- gsub("^#", "", cadd_annotations_2_header)
cadd_annotations_2_header <- strsplit(cadd_annotations_2_header, split="\t")[[1]]

cadd_annotations_2 <- read.delim(cadd_annotations_2_saveFile, header=F, sep="\t", stringsAsFactors=F, comment.char="#")
colnames(cadd_annotations_2) <- cadd_annotations_2_header

#
cadd_annotations_2 <- cadd_annotations_2[,c("Location", "GerpRS", "CADD_RAW", "CADD_PHRED")]

cadd_annotations_2_locations <- as.data.frame(matrix(nrow=nrow(cadd_annotations_2), ncol=2, data=unlist(strsplit(as.character(cadd_annotations_2$Location), split=":", fixed=T)), byrow=T))
colnames(cadd_annotations_2_locations) <- c("chr", "start")
withDash <- grep("-", cadd_annotations_2_locations$start)
cadd_annotations_2_locations[-withDash, "start"] <- paste(cadd_annotations_2_locations[-withDash, "start"], cadd_annotations_2_locations[-withDash, "start"], sep="-")
corrected <- as.data.frame(matrix(nrow=nrow(cadd_annotations_2_locations), ncol=2, data=unlist(strsplit(cadd_annotations_2_locations[,2], split="-")), byrow=T))
cadd_annotations_2_locations <- cbind(cadd_annotations_2_locations$chr, corrected)
colnames(cadd_annotations_2_locations) <- c("chr", "start", "end")

cadd_annotations_2 <- cbind(cadd_annotations_2_locations, cadd_annotations_2[,2:4])
cadd_annotations_2_gr <- makeGRangesFromDataFrame(cadd_annotations_2)










forbiddenExprs <- colnames(pval_table_summed_1)[c(grep("H[3|4]K", colnames(pval_table_summed_1)), grep("INPUT", colnames(pval_table_summed_1)), grep("^POL", colnames(pval_table_summed_1)))]
pval_table_summed_1 <- pval_table_summed_1[,!colnames(pval_table_summed_1)%in%forbiddenExprs]
depth_table_summed_1 <- depth_table_summed_1[,!colnames(depth_table_summed_1)%in%forbiddenExprs]

toRemove <- grep("chrX", rownames(pval_table_summed_1))
if(length(toRemove)>0) {
  pval_table_summed_1 <- pval_table_summed_1[-toRemove,]
  depth_table_summed_1 <- depth_table_summed_1[-toRemove,]
}


forbiddenExprs <- colnames(pval_table_summed_2)[c(grep("H[3|4]K", colnames(pval_table_summed_2)), grep("INPUT", colnames(pval_table_summed_2)), grep("^POL", colnames(pval_table_summed_2)))]
pval_table_summed_2 <- pval_table_summed_2[,!colnames(pval_table_summed_2)%in%forbiddenExprs]
depth_table_summed_2 <- depth_table_summed_2[,!colnames(depth_table_summed_2)%in%forbiddenExprs]


toRemove <- grep("chrX", rownames(pval_table_summed_2))
if(length(toRemove)>0) {
  pval_table_summed_2 <- pval_table_summed_2[-toRemove,]
  depth_table_summed_2 <- depth_table_summed_2[-toRemove,]
}




regionsMat_d1 <- matrix(nrow=nrow(pval_table_summed_1), ncol=3, data=unlist(strsplit(rownames(pval_table_summed_1), split=":|-")), byrow=T)
regionsMat_d2 <- matrix(nrow=nrow(pval_table_summed_2), ncol=3, data=unlist(strsplit(rownames(pval_table_summed_2), split=":|-")), byrow=T)

regionsMat_d1 <- as.data.frame(regionsMat_d1)
regionsMat_d1[,1] <- as.character(regionsMat_d1[,1])
regionsMat_d1[,2] <- as.numeric(as.character(regionsMat_d1[,2]))
regionsMat_d1[,3] <- as.numeric(as.character(regionsMat_d1[,3]))
colnames(regionsMat_d1) <- c("chr", "start", "end")
regionsMat_d1_gr <- makeGRangesFromDataFrame(regionsMat_d1, ignore.strand=T, keep.extra.columns=T)

regionsMat_d2 <- as.data.frame(regionsMat_d2)
regionsMat_d2[,1] <- as.character(regionsMat_d2[,1])
regionsMat_d2[,2] <- as.numeric(as.character(regionsMat_d2[,2]))
regionsMat_d2[,3] <- as.numeric(as.character(regionsMat_d2[,3]))
colnames(regionsMat_d2) <- c("chr", "start", "end")
regionsMat_d2_gr <- makeGRangesFromDataFrame(regionsMat_d2, ignore.strand=T, keep.extra.columns=T)






##################################################################
#Find appropriate overlaps for all relevant data, and get enrichments.
##################################################################




conservation_overlaps_d1 <- as.data.frame(findOverlaps(regionsMat_d1_gr, cadd_annotations_1_gr))
conservation_overlaps_d2 <- as.data.frame(findOverlaps(regionsMat_d2_gr, cadd_annotations_2_gr))



enrichmnt_Table_11_combined <- getEnrichmentTable_combined(pval_table_summed_1, depth_table_summed_1, conservation_overlaps_d1, cadd_annotations_1, pval_table_summed_2, depth_table_summed_2, conservation_overlaps_d2, cadd_annotations_2, suffDepth=11)




#Save a table of those.

saveFile <- paste(outDir, "Enrichment_Table_combinedDonors_depth11.txt", sep="")
write.table(enrichmnt_Table_11_combined, saveFile, row.names=F, col.names=T, sep="\t", quote=F)



#Save a table of each of those.

saveBase <- paste(outDir, "Enrichment_Barplot_combined_depth11", sep="")
plot_enrichments(enrichmnt_Table_11_combined, saveBase)





#
