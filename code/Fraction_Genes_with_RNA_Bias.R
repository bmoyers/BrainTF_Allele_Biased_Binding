#Fraction_Genes_with_RNA_Bias.R



################################################################################
#This script is used for minor analysis to determine the fraction of allele-biased
#     cases of RNA seq which fall within accepted gene models, and the fraction
#     of gene models which show allele bias.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges, matrixStats, and TxDb.Hsapiens.UCSC.hg38.knownGene packages
#     are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# rna_d1: A table compiling p-value calculations for RNA in donor 1, produced
#     by the script CompilingHaplotype_RNA_tables_summedAcrossTissues.R.
# rna_d2: As rna_d1, but for donor 2.
#
################################################################################



#cd ~/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/
#module load cluster/R/4.1.0

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(matrixStats)
library(GenomicRanges)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")


################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
rna_d1 <- args[2]
rna_d2 <- args[3]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#rna_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_RNA_binom_vg.txt"
#rna_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_RNA_binom_vg.txt"



################################################################################
#Get a genomicRanges of all gene bodies.
################################################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

geneBodies <- GenomicFeatures::genes(txdb)




################################################################################
#Load in the RNA-pval tables from each donor.
#Determine, for each donor, the fraction of gene bodies overlapping with
#an Allele-biased expression at P<=0.001.
################################################################################


rna_d1 <- read.table(rna_d1, header=T, sep="\t", stringsAsFactors=F)
rna_d1 <- rna_d1[rna_d1$pvals_summed<=0.001,]
rna_d1_regions <- as.data.frame(matrix(nrow=nrow(rna_d1), ncol=3, data=unlist(strsplit(rownames(rna_d1), split=":|-")), byrow=T))
colnames(rna_d1_regions) <- c("chr", "start", "end")
rna_d1_regions_gr <- makeGRangesFromDataFrame(rna_d1_regions)

theIntersect <- as.data.frame(findOverlaps(rna_d1_regions_gr, geneBodies))
fraction_d1_ase <- length(unique(theIntersect[,2])) / length(geneBodies)

rna_d2 <- read.table(rna_d2, header=T, sep="\t", stringsAsFactors=F)
rna_d2 <- rna_d2[rna_d2$pvals_summed<=0.001,]
rna_d2_regions <- as.data.frame(matrix(nrow=nrow(rna_d2), ncol=3, data=unlist(strsplit(rownames(rna_d2), split=":|-")), byrow=T))
colnames(rna_d2_regions) <- c("chr", "start", "end")
rna_d2_regions_gr <- makeGRangesFromDataFrame(rna_d2_regions)

theIntersect <- as.data.frame(findOverlaps(rna_d2_regions_gr, geneBodies))
fraction_d2_ase <- length(unique(theIntersect[,2])) / length(geneBodies)











#
