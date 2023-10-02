#!/usr/bin/env R
#Supplemental_Table_5.R

################################################################################
#This script is used to produce Supplemental Table 4.
#
#Run under R 4.1.0, but can be used under other versions so long as the packages
#     GenomicRanges, matrixStats, and ggplot are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# gtex_dir: Directory containing the tissue eQTLs from GTEx, downloaded from
#     https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
#     on 30 June 2023.
# d1_sequences: path to the fasta file for all variant sequences in donor 1,
#     produced by the script Building_fasta_variant_sequences.sh
# pval_table_1: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# d2_sequences: As d1_sequences, but for donor 2
# pval_table_2: As pval_table_1, but for donor 2
#
################################################################################

################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(matrixStats)
library(GenomicRanges)

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
gtex_dir <- args[2]
d1_sequences <- args[3]
pval_table_1 <- args[4]
d2_sequences <- args[5]
pval_table_2 <- args[6]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#gtex_dir <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/GTEx_Analysis_v8_eQTL/"
#d1_sequences <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#d2_sequences <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"



################################################################################
#Load in GTEx data.  downloaded from:
#https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
#on 30 June 2023.
################################################################################

#brainFiles <- list.files(gtex_dir, full.names=T, pattern="Brain")
brainFiles <- list.files(gtex_dir, full.names=T)
brainFiles_pairs <- brainFiles[grep("signif_variant_gene_pairs", brainFiles)]
brainFiles_eGenes <- brainFiles[-grep("signif_variant_gene_pairs", brainFiles)]


fullGTEX <- c()
for (i in 1:length(brainFiles_pairs)) {
  print(paste(i, length(brainFiles_pairs)))
  thesePairs <- read.table(gzfile(brainFiles_pairs[i]), header=T, sep="\t", stringsAsFactors=F)
  theseEgenes <- read.table(gzfile(brainFiles_eGenes[i]), header=T, sep="\t", stringsAsFactors=F)

  theseEgenes_match <- theseEgenes[match(thesePairs$gene_id, theseEgenes$gene_id),c("gene_name", "gene_chr", "gene_start", "gene_end", "strand")]
  thesePairs <- cbind(thesePairs, theseEgenes_match)
  thisTissue <- strsplit(brainFiles_pairs[i], split="/")[[1]]
  thisTissue <- thisTissue[length(thisTissue)]
  thisTissue <- strsplit(thisTissue, ".v8.", fixed=T)[[1]][1]
  thesePairs$Tissue <- rep(thisTissue, nrow(thesePairs))
  these_variants <- unique(thesePairs$variant_id)
  these_variants <- these_variants[!these_variants%in%fullGTEX]
  fullGTEX <- c(fullGTEX, these_variants)
}

fullGTEX <- fullGTEX[order(fullGTEX)]
fullGTEX_regions <- as.data.frame(matrix(nrow=length(fullGTEX), ncol=5, data=unlist(strsplit(fullGTEX, split="_")), byrow=T))
fullGTEX_regions <- fullGTEX_regions[,c(1,2,2,3,4,5)]
colnames(fullGTEX_regions) <- c("chr", "start", "end", "a1", "a2", "etc")
fullGTEX_regions_gr <- makeGRangesFromDataFrame(fullGTEX_regions)

fullGTEX_regions$altNames <- paste(fullGTEX_regions[,1], fullGTEX_regions[,2], fullGTEX_regions[,4], fullGTEX_regions[,5], sep="_")
fullGTEX_regions$hom1 <- paste(fullGTEX_regions[,1], fullGTEX_regions[,2], fullGTEX_regions[,4], sep="_")
fullGTEX_regions$hom2 <- paste(fullGTEX_regions[,1], fullGTEX_regions[,2], fullGTEX_regions[,5], sep="_")


################################################################################
################################################################################
#Load in the variants for each donor.
#Determine, for each eQTL in GTEx, whether or not the variant
#is present in donor 1 (and in what state), and donor 2 (and in what state)
################################################################################
################################################################################

################################################################################
################################################################################
#Donor 1
################################################################################
################################################################################


d1_sequences <- readLines(d1_sequences)
d1_sequences <- d1_sequences[grep("^>",d1_sequences)]
d1_sequences <- gsub(">", "", d1_sequences)

d1_sequences_mat <- d1_sequences[grep(":M:", d1_sequences)]
d1_sequences_pat <- d1_sequences[grep(":P:", d1_sequences)]

d1_sequences_mat <- matrix(nrow=length(d1_sequences_mat), ncol=2, data=unlist(strsplit(d1_sequences_mat, split=":M:", fixed=T)), byrow=T)
d1_sequences_pat <- matrix(nrow=length(d1_sequences_pat), ncol=2, data=unlist(strsplit(d1_sequences_pat, split=":P:", fixed=T)), byrow=T)


d1_regions_w_vars <- matrix(nrow=nrow(d1_sequences_mat), ncol=3, data=unlist(strsplit(d1_sequences_mat[,1], split=":|-")), byrow=T)
d1_regions_w_vars <- as.data.frame(d1_regions_w_vars)
colnames(d1_regions_w_vars) <- c("chr", "start", "end")
d1_regions_w_vars$region <- d1_sequences_mat[,1]
d1_regions_w_vars$matVar <- d1_sequences_mat[,2]
d1_regions_w_vars$patVar <- d1_sequences_pat[,2]

d1_regions_w_vars$start <- as.numeric(as.character(d1_regions_w_vars$start))
d1_regions_w_vars$end <- as.numeric(as.character(d1_regions_w_vars$end))

d1_regions_w_vars$start <- d1_regions_w_vars$start+100
d1_regions_w_vars$end <- d1_regions_w_vars$end-100
d1_regions_w_vars_gr <- makeGRangesFromDataFrame(d1_regions_w_vars)

theIntersect <- as.data.frame(findOverlaps(d1_regions_w_vars_gr, fullGTEX_regions_gr))

d1_regions_w_vars <- d1_regions_w_vars[unique(theIntersect[,1]),]



d1_state <- rep("HomozygousReference", length(fullGTEX))
d1_significant <- rep("Nonsignificant", length(fullGTEX))


################################################################################
#Determine the state of each variant-- Het, Homozygous Reference,
#or Homozygous Alternate.
################################################################################

allVars_mat <- unlist(strsplit(d1_regions_w_vars$matVar, split=";"))

allVars_mat_het <- allVars_mat[grep("HET", allVars_mat)]
matVars_het_length <- length(allVars_mat_het)

d1_vars_alternateFormat_Hets <- matrix(nrow=matVars_het_length, ncol=5, data=NA)
currentRow <- 1
for (i in 1:nrow(d1_regions_w_vars)) {
  if(i%%1000==0) {print(paste(i, nrow(d1_regions_w_vars)))}
  thisSet_mat <- unlist(strsplit(d1_regions_w_vars[i,"matVar"], split=";"))
  thisSet_mat <- as.data.frame(matrix(nrow=length(thisSet_mat), ncol=3, data=unlist(strsplit(thisSet_mat, split=",")), byrow=T))
  thisSet_pat <- unlist(strsplit(d1_regions_w_vars[i,"patVar"], split=";"))
  thisSet_pat <- as.data.frame(matrix(nrow=length(thisSet_pat), ncol=3, data=unlist(strsplit(thisSet_pat, split=",")), byrow=T))
  thisSet_mat <- thisSet_mat[thisSet_mat[,3]=="HET",]
  thisSet_pat <- thisSet_pat[thisSet_pat[,3]=="HET",]
  if(nrow(thisSet_mat)>0) {
    for (j in 1:nrow(thisSet_mat)) {
      thisLine <- c(d1_regions_w_vars[i,"chr"], as.character(thisSet_mat[j,1:2]), as.character(thisSet_pat[j,2:3]))
      d1_vars_alternateFormat_Hets[currentRow,] <- thisLine
      currentRow <- currentRow+1
    }
  }

}

d1_vars_alternateFormat_Hets <- as.data.frame(d1_vars_alternateFormat_Hets)
colnames(d1_vars_alternateFormat_Hets) <- c("chr", "loc", "mat", "pat", "type")
d1_vars_alternateFormat_Hets$name1 <- paste(d1_vars_alternateFormat_Hets[,1], d1_vars_alternateFormat_Hets[,2], d1_vars_alternateFormat_Hets[,3], d1_vars_alternateFormat_Hets[,4], sep="_")
d1_vars_alternateFormat_Hets$name2 <- paste(d1_vars_alternateFormat_Hets[,1], d1_vars_alternateFormat_Hets[,2], d1_vars_alternateFormat_Hets[,4], d1_vars_alternateFormat_Hets[,3], sep="_")


d1_state[which(fullGTEX_regions$altNames%in%d1_vars_alternateFormat_Hets$name1)] <- "Heterozygous"
d1_state[which(fullGTEX_regions$altNames%in%d1_vars_alternateFormat_Hets$name2)] <- "Heterozygous"



allVars_mat_hom <- allVars_mat[grep("HOM", allVars_mat)]
matVars_hom_length <- length(allVars_mat_hom)

d1_vars_alternateFormat_Homs <- matrix(nrow=matVars_hom_length, ncol=5, data=NA)
currentRow <- 1
for (i in 1:nrow(d1_regions_w_vars)) {
  if(i%%1000==0) {print(paste(i, nrow(d1_regions_w_vars)))}
  thisSet_mat <- unlist(strsplit(d1_regions_w_vars[i,"matVar"], split=";"))
  thisSet_mat <- as.data.frame(matrix(nrow=length(thisSet_mat), ncol=3, data=unlist(strsplit(thisSet_mat, split=",")), byrow=T))
  thisSet_pat <- unlist(strsplit(d1_regions_w_vars[i,"patVar"], split=";"))
  thisSet_pat <- as.data.frame(matrix(nrow=length(thisSet_pat), ncol=3, data=unlist(strsplit(thisSet_pat, split=",")), byrow=T))
  thisSet_mat <- thisSet_mat[thisSet_mat[,3]=="HOM",]
  thisSet_pat <- thisSet_pat[thisSet_pat[,3]=="HOM",]
  if(nrow(thisSet_mat)>0) {
    for (j in 1:nrow(thisSet_mat)) {
      thisLine <- c(d1_regions_w_vars[i,"chr"], as.character(thisSet_mat[j,1:2]), as.character(thisSet_pat[j,2:3]))
      d1_vars_alternateFormat_Homs[currentRow,] <- thisLine
      currentRow <- currentRow+1
    }
  }

}

d1_vars_alternateFormat_Homs <- as.data.frame(d1_vars_alternateFormat_Homs)
colnames(d1_vars_alternateFormat_Homs) <- c("chr", "loc", "mat", "pat", "type")
d1_vars_alternateFormat_Homs$name <- paste(d1_vars_alternateFormat_Homs[,1], d1_vars_alternateFormat_Homs[,2], d1_vars_alternateFormat_Homs[,3], sep="_")

a1 <- which(fullGTEX_regions$hom1%in%d1_vars_alternateFormat_Homs$name)
a2 <- which(fullGTEX_regions$hom2%in%d1_vars_alternateFormat_Homs$name)
a <- c(a1, a2)
d1_state[a] <- "HomozygousAlternate"



fullGTEX_regions$d1_state <- d1_state




################################################################################
#Determine which variants overlap with a region that has allele-biased
#binding.  NOTE that some homozygous variants might have this because there
#is a closely-assocaited SNP which is het and has ABB.
#So for the sake of simplicity, let's only do this for het cases.
################################################################################


pval_table_1 <- read.table(pval_table_1, header=T, sep="\t", stringsAsFactors=F)

significant_input_d1 <- pval_table_1[pval_table_1[,"INPUT"]<=0.05,]

toRemove <- grep("chrX", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
}

forbiddenExprs <- colnames(pval_table_1)[c(grep("H[3|4]K", colnames(pval_table_1)), grep("INPUT", colnames(pval_table_1)), grep("POL2", colnames(pval_table_1)))]
pval_table_1 <- pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs]

pval_table_1 <- pval_table_1[rowMins(as.matrix(pval_table_1))<=0.001,]
pval_table_1_regions <- as.data.frame(matrix(nrow=nrow(pval_table_1), ncol=3, data=unlist(strsplit(rownames(pval_table_1), split=":|-")), byrow=T))
colnames(pval_table_1_regions) <- c("chr", "start", "end")
pval_table_1_regions <- makeGRangesFromDataFrame(pval_table_1_regions)


theIntersect <- as.data.frame(findOverlaps(fullGTEX_regions_gr, pval_table_1_regions))
areSig <- unique(theIntersect[,1])
areHet <- which(d1_state=="Heterozygous")
areSig <- areSig[areSig%in%areHet]
d1_significant[areSig] <- "Significant"


fullGTEX_regions$d1_tf_significance <- d1_significant


#




################################################################################
################################################################################
#Donor 2
################################################################################
################################################################################


d2_sequences <- readLines(d2_sequences)
d2_sequences <- d2_sequences[grep("^>",d2_sequences)]
d2_sequences <- gsub(">", "", d2_sequences)

d2_sequences_mat <- d2_sequences[grep(":M:", d2_sequences)]
d2_sequences_pat <- d2_sequences[grep(":P:", d2_sequences)]

d2_sequences_mat <- matrix(nrow=length(d2_sequences_mat), ncol=2, data=unlist(strsplit(d2_sequences_mat, split=":M:", fixed=T)), byrow=T)
d2_sequences_pat <- matrix(nrow=length(d2_sequences_pat), ncol=2, data=unlist(strsplit(d2_sequences_pat, split=":P:", fixed=T)), byrow=T)


d2_regions_w_vars <- matrix(nrow=nrow(d2_sequences_mat), ncol=3, data=unlist(strsplit(d2_sequences_mat[,1], split=":|-")), byrow=T)
d2_regions_w_vars <- as.data.frame(d2_regions_w_vars)
colnames(d2_regions_w_vars) <- c("chr", "start", "end")
d2_regions_w_vars$region <- d2_sequences_mat[,1]
d2_regions_w_vars$matVar <- d2_sequences_mat[,2]
d2_regions_w_vars$patVar <- d2_sequences_pat[,2]

d2_regions_w_vars$start <- as.numeric(as.character(d2_regions_w_vars$start))
d2_regions_w_vars$end <- as.numeric(as.character(d2_regions_w_vars$end))

d2_regions_w_vars$start <- d2_regions_w_vars$start+100
d2_regions_w_vars$end <- d2_regions_w_vars$end-100
d2_regions_w_vars_gr <- makeGRangesFromDataFrame(d2_regions_w_vars)

theIntersect <- as.data.frame(findOverlaps(d2_regions_w_vars_gr, fullGTEX_regions_gr))

d2_regions_w_vars <- d2_regions_w_vars[unique(theIntersect[,1]),]



d2_state <- rep("HomozygousReference", length(fullGTEX))
d2_significant <- rep("Nonsignificant", length(fullGTEX))


################################################################################
#Determine the state of each variant-- Het, Homozygous Reference,
#or Homozygous Alternate.
################################################################################

allVars_mat <- unlist(strsplit(d2_regions_w_vars$matVar, split=";"))

allVars_mat_het <- allVars_mat[grep("HET", allVars_mat)]
matVars_het_length <- length(allVars_mat_het)

d2_vars_alternateFormat_Hets <- matrix(nrow=matVars_het_length, ncol=5, data=NA)
currentRow <- 1
for (i in 1:nrow(d2_regions_w_vars)) {
  if(i%%1000==0) {print(paste(i, nrow(d2_regions_w_vars)))}
  thisSet_mat <- unlist(strsplit(d2_regions_w_vars[i,"matVar"], split=";"))
  thisSet_mat <- as.data.frame(matrix(nrow=length(thisSet_mat), ncol=3, data=unlist(strsplit(thisSet_mat, split=",")), byrow=T))
  thisSet_pat <- unlist(strsplit(d2_regions_w_vars[i,"patVar"], split=";"))
  thisSet_pat <- as.data.frame(matrix(nrow=length(thisSet_pat), ncol=3, data=unlist(strsplit(thisSet_pat, split=",")), byrow=T))
  thisSet_mat <- thisSet_mat[thisSet_mat[,3]=="HET",]
  thisSet_pat <- thisSet_pat[thisSet_pat[,3]=="HET",]
  if(nrow(thisSet_mat)>0) {
    for (j in 1:nrow(thisSet_mat)) {
      thisLine <- c(d2_regions_w_vars[i,"chr"], as.character(thisSet_mat[j,1:2]), as.character(thisSet_pat[j,2:3]))
      d2_vars_alternateFormat_Hets[currentRow,] <- thisLine
      currentRow <- currentRow+1
    }
  }

}

d2_vars_alternateFormat_Hets <- as.data.frame(d2_vars_alternateFormat_Hets)
colnames(d2_vars_alternateFormat_Hets) <- c("chr", "loc", "mat", "pat", "type")
d2_vars_alternateFormat_Hets$name1 <- paste(d2_vars_alternateFormat_Hets[,1], d2_vars_alternateFormat_Hets[,2], d2_vars_alternateFormat_Hets[,3], d2_vars_alternateFormat_Hets[,4], sep="_")
d2_vars_alternateFormat_Hets$name2 <- paste(d2_vars_alternateFormat_Hets[,1], d2_vars_alternateFormat_Hets[,2], d2_vars_alternateFormat_Hets[,4], d2_vars_alternateFormat_Hets[,3], sep="_")


d2_state[which(fullGTEX_regions$altNames%in%d2_vars_alternateFormat_Hets$name1)] <- "Heterozygous"
d2_state[which(fullGTEX_regions$altNames%in%d2_vars_alternateFormat_Hets$name2)] <- "Heterozygous"



allVars_mat_hom <- allVars_mat[grep("HOM", allVars_mat)]
matVars_hom_length <- length(allVars_mat_hom)

d2_vars_alternateFormat_Homs <- matrix(nrow=matVars_hom_length, ncol=5, data=NA)
currentRow <- 1
for (i in 1:nrow(d2_regions_w_vars)) {
  if(i%%1000==0) {print(paste(i, nrow(d2_regions_w_vars)))}
  thisSet_mat <- unlist(strsplit(d2_regions_w_vars[i,"matVar"], split=";"))
  thisSet_mat <- as.data.frame(matrix(nrow=length(thisSet_mat), ncol=3, data=unlist(strsplit(thisSet_mat, split=",")), byrow=T))
  thisSet_pat <- unlist(strsplit(d2_regions_w_vars[i,"patVar"], split=";"))
  thisSet_pat <- as.data.frame(matrix(nrow=length(thisSet_pat), ncol=3, data=unlist(strsplit(thisSet_pat, split=",")), byrow=T))
  thisSet_mat <- thisSet_mat[thisSet_mat[,3]=="HOM",]
  thisSet_pat <- thisSet_pat[thisSet_pat[,3]=="HOM",]
  if(nrow(thisSet_mat)>0) {
    for (j in 1:nrow(thisSet_mat)) {
      thisLine <- c(d2_regions_w_vars[i,"chr"], as.character(thisSet_mat[j,1:2]), as.character(thisSet_pat[j,2:3]))
      d2_vars_alternateFormat_Homs[currentRow,] <- thisLine
      currentRow <- currentRow+1
    }
  }

}

d2_vars_alternateFormat_Homs <- as.data.frame(d2_vars_alternateFormat_Homs)
colnames(d2_vars_alternateFormat_Homs) <- c("chr", "loc", "mat", "pat", "type")
d2_vars_alternateFormat_Homs$name <- paste(d2_vars_alternateFormat_Homs[,1], d2_vars_alternateFormat_Homs[,2], d2_vars_alternateFormat_Homs[,3], sep="_")

a1 <- which(fullGTEX_regions$hom1%in%d2_vars_alternateFormat_Homs$name)
a2 <- which(fullGTEX_regions$hom2%in%d2_vars_alternateFormat_Homs$name)
a <- c(a1, a2)
d2_state[a] <- "HomozygousAlternate"



fullGTEX_regions$d2_state <- d2_state




################################################################################
#Determine which variants overlap with a region that has allele-biased
#binding.  NOTE that some homozygous variants might have this because there
#is a closely-assocaited SNP which is het and has ABB.
#So for the sake of simplicity, let's only do this for het cases.
################################################################################


pval_table_2 <- read.table(pval_table_2, header=T, sep="\t", stringsAsFactors=F)

significant_input_d2 <- pval_table_2[pval_table_2[,"INPUT"]<=0.05,]

toRemove <- grep("chrX", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
}

forbiddenExprs <- colnames(pval_table_2)[c(grep("H[3|4]K", colnames(pval_table_2)), grep("INPUT", colnames(pval_table_2)), grep("POL2", colnames(pval_table_2)))]
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs]

pval_table_2 <- pval_table_2[rowMins(as.matrix(pval_table_2))<=0.001,]
pval_table_2_regions <- as.data.frame(matrix(nrow=nrow(pval_table_2), ncol=3, data=unlist(strsplit(rownames(pval_table_2), split=":|-")), byrow=T))
colnames(pval_table_2_regions) <- c("chr", "start", "end")
pval_table_2_regions <- makeGRangesFromDataFrame(pval_table_2_regions)


theIntersect <- as.data.frame(findOverlaps(fullGTEX_regions_gr, pval_table_2_regions))
areSig <- unique(theIntersect[,1])
areHet <- which(d2_state=="Heterozygous")
areSig <- areSig[areSig%in%areHet]
d2_significant[areSig] <- "Significant"


fullGTEX_regions$d2_tf_significance <- d2_significant


#



################################################################################
################################################################################
################################################################################
#Restrict to cases that are present in at least one of the two donors.
#Save the Table
################################################################################
################################################################################
################################################################################

set1 <- which(fullGTEX_regions$d1_state!="HomozygousReference")
set2 <- which(fullGTEX_regions$d2_state!="HomozygousReference")
setBoth <- unique(c(set1, set2))
setBoth <- setBoth[order(setBoth)]

fullGTEX_regions_inDonors <- fullGTEX_regions[setBoth,]


set1 <- which(fullGTEX_regions$d1_tf_significance=="Significant")
set2 <- which(fullGTEX_regions$d2_tf_significance=="Significant")
setBoth <- unique(c(set1, set2))
setBoth <- setBoth[order(setBoth)]

fullGTEX_regions_inDonors <- fullGTEX_regions[setBoth,]


saveFile <- paste(outDir, "Supplemental_Table_5.txt", sep="")
write.table(fullGTEX_regions_inDonors, saveFile, row.names=F, col.names=T, sep="\t", quote=F)

#
