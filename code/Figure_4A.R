#!/usr/bin/env R
#Figure_4A.R


################################################################################
#This script is used to produce Figure 4A.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges, ggplot2, and matrixStats packages are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# gtex_dir: Directory containing the tissue eQTLs from GTEx, downloaded from
#     https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
#     on 30 June 2023.
# rna_tf_sig_saveFile_d1: A table compiling cases of significant allele-biased RNA
#     expression for donor 1, produced by the script
# d1_sequences: A fasta file of all variant sequences for donor 1, produced
#     by the script Building_fasta_variant_sequences.sh .
# pval_table_1: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# rna_tf_sig_saveFile_d2: As rna_tf_sig_saveFile_d1, but for donor 2.
# d2_sequences: As d1_sequences, but for donor 2.
# pval_table_2: As pval_table_1, but for donor 2.
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

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

makeBarplotGraphDF <- function(RNA_table) {
  #RNA_table <- donor1_relevant_RNA
  ref_df <- cbind(RNA_table$rna_refReads, RNA_table$rna_region, rep("Reference", nrow(RNA_table)))
  colnames(ref_df) <- c("Reads", "Region", "type")
  alt_df <- cbind(RNA_table$rna_altReads, RNA_table$rna_region, rep("Alternate", nrow(RNA_table)))
  colnames(ref_df) <- c("Reads", "Region", "type")

  finalDF <- rbind(ref_df, alt_df)
  finalDF <- as.data.frame(finalDF)
  finalDF[,1] <- as.numeric(as.character(finalDF[,1]))
  #finalDF[,2] <- factor(finalDF[,2], levels=unique(finalDF[,2]))
  finalDF[,2] <- as.character(finalDF[,2])
  finalDF[,3] <- factor(finalDF[,3], levels=unique(finalDF[,3]))
  return(finalDF)
}


################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
gtex_dir <- args[2]
rna_tf_sig_saveFile_d1 <- args[3]
d1_sequences <- args[4]
pval_table_1 <- args[5]
rna_tf_sig_saveFile_d2 <- args[6]
d2_sequences <- args[7]
pval_table_2 <- args[8]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#gtex_dir <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/GTEx_Analysis_v8_eQTL/"
#rna_tf_sig_saveFile_d1 <- paste(outDir, "rna_sig_d1_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
#d1_sequences <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#rna_tf_sig_saveFile_d2 <- paste(outDir, "rna_sig_d2_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
#d2_sequences <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"



################################################################################
#Load in GTEx data.  downloaded from:
#https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
#on 30 June 2023.
################################################################################

brainFiles <- list.files(gtex_dir, full.names=T, pattern="Brain")
#brainFiles <- list.files(gtex_dir, full.names=T)
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
  fullGTEX <- rbind(fullGTEX, thesePairs)
}

fullGTEX$alt_variant_id <- fullGTEX$variant_id
fullGTEX$alt_variant_id <- gsub("_b38", "", fullGTEX$alt_variant_id)



GTEX_regions <- matrix(nrow=nrow(fullGTEX), ncol=5, data=unlist(strsplit(fullGTEX[,1], split="_")), byrow=T)
GTEX_regions[,3] <- GTEX_regions[,2]
GTEX_regions <- GTEX_regions[,1:3]
colnames(GTEX_regions) <- c("chr", "start", "end")
GTEX_regions <- makeGRangesFromDataFrame(as.data.frame(GTEX_regions))



################################################################################
################################################################################
#Donor1
################################################################################
################################################################################


################################################################################
#Load in linked RNA-TF bias tables.  Use TF information there to restrict
#GTEx table to variants within a given donor.
################################################################################


rna_tf_sig_wPhase_d1 <- read.table(rna_tf_sig_saveFile_d1, header=T, sep="\t", stringsAsFactors=F)
relevant_tf_regions_d1 <- rna_tf_sig_wPhase_d1[,4:6]
colnames(relevant_tf_regions_d1) <- c("chr", "start", "end")
relevant_tf_regions_d1 <- makeGRangesFromDataFrame(relevant_tf_regions_d1)
theIntersect <- as.data.frame(findOverlaps(GTEX_regions, relevant_tf_regions_d1))

d1_fullGTEX <- fullGTEX[unique(theIntersect[,1]),]
rna_tf_sig_wPhase_d1 <- rna_tf_sig_wPhase_d1[unique(theIntersect[,2]),]

################################################################################
#Restrict both gtex variants and our rna variants to cases where the
#RNA variation is within the appropriate gene.
################################################################################

d1_fullGTEX_rna_regions <- d1_fullGTEX[,c("gene_chr", "gene_start", "gene_end")]
colnames(d1_fullGTEX_rna_regions) <- c("chr", "start", "end")
d1_fullGTEX_rna_regions <- makeGRangesFromDataFrame(d1_fullGTEX_rna_regions)

relevant_rna_regions_d1 <- rna_tf_sig_wPhase_d1[,1:3]
colnames(relevant_rna_regions_d1) <- c("chr", "start", "end")
relevant_rna_regions_d1 <- makeGRangesFromDataFrame(relevant_rna_regions_d1)

theIntersect <- as.data.frame(findOverlaps(d1_fullGTEX_rna_regions, relevant_rna_regions_d1))

d1_fullGTEX <- d1_fullGTEX[unique(theIntersect[,1]),]
rna_tf_sig_wPhase_d1 <- rna_tf_sig_wPhase_d1[unique(theIntersect[,2]),]


################################################################################
#For ease in the next section, restrict the gtex data and
#add a tf region column to the rna_tf_sig_wPhase table
################################################################################

d1_fullGTEX <- d1_fullGTEX[,c("alt_variant_id", "gene_id", "tss_distance", "maf", "pval_nominal", "slope", "gene_name", "gene_chr", "gene_start", "gene_end", "strand", "Tissue")]

rna_tf_sig_wPhase_d1$tf_region <- paste(rna_tf_sig_wPhase_d1[,"chr_tf"], ":", rna_tf_sig_wPhase_d1[,"start_tf"], "-", rna_tf_sig_wPhase_d1[,"end_tf"], sep="")
rna_tf_sig_wPhase_d1$rna_region <- paste(rna_tf_sig_wPhase_d1[,"chr"], ":", rna_tf_sig_wPhase_d1[,"start"], "-", rna_tf_sig_wPhase_d1[,"end"], sep="")

################################################################################
#We need to restrict the table further based on actual
#variants that we observe.  To do this, I need to read in the
#names of all variant regions and parse out, for each GTEx entry:
#1) do we have a variant at the relevant location
#2) if so, do we have the appropriate reference and alternate allele?
#3) if so, correctly sort the read counts from mat/pat to ref/alt for
#   a proper comparison to GTEx predictions.
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

d1_relevant_regions <- paste(rna_tf_sig_wPhase_d1[,4], ":", rna_tf_sig_wPhase_d1[,5], "-", rna_tf_sig_wPhase_d1[,6], sep="")
d1_regions_w_vars <- d1_regions_w_vars[d1_regions_w_vars$region%in%d1_relevant_regions,]
d1_regions_w_vars_gr <- makeGRangesFromDataFrame(d1_regions_w_vars)



d1_fullGTEX_inOurSet <- c()

for (i in 1:nrow(d1_fullGTEX)) {
  if(i%%100==0) {print(paste(i, nrow(d1_fullGTEX), nrow(d1_fullGTEX_inOurSet)))}
  thisVar <- strsplit(d1_fullGTEX[i,"alt_variant_id"], split="_")[[1]]
  thisVar_gr <- data.frame(chr=thisVar[1], start=as.numeric(thisVar[2]), end=as.numeric(thisVar[2]))
  thisVar_gr <- makeGRangesFromDataFrame(thisVar_gr)

  theIntersect <- as.data.frame(findOverlaps(d1_regions_w_vars_gr, thisVar_gr))
  if(nrow(theIntersect)>0) {
    tf_region <- d1_regions_w_vars[theIntersect[1,1],"region"]

    matVar <- d1_regions_w_vars[theIntersect[1,1],"matVar"]
    patVar <- d1_regions_w_vars[theIntersect[1,1],"patVar"]

    matVar <- strsplit(matVar, split=";")[[1]]
    patVar <- strsplit(patVar, split=";")[[1]]

    matVar <- matrix(nrow=length(matVar), ncol=3, data=unlist(strsplit(matVar, split=",")), byrow=T)
    patVar <- matrix(nrow=length(patVar), ncol=3, data=unlist(strsplit(patVar, split=",")), byrow=T)

    matVar <- as.data.frame(matVar)
    colnames(matVar) <- c("loc", "nuc", "type")
    patVar <- as.data.frame(patVar)
    colnames(patVar) <- c("loc", "nuc", "type")

    matVar <- matVar[matVar[,3]=="HET",]
    patVar <- patVar[patVar[,3]=="HET",]

    matVar <- matVar[matVar[,"loc"]==thisVar[2],]
    patVar <- patVar[patVar[,"loc"]==thisVar[2],]

    if(nrow(matVar)>1 || nrow(patVar)>1) {
      print(i)
      break()
    }

    if(nrow(matVar)==1 && nrow(patVar)==1) {
      thisLine_GTEX <- d1_fullGTEX[i,]
      thisSet <- rna_tf_sig_wPhase_d1[rna_tf_sig_wPhase_d1$tf_region==tf_region,]

      thisLine_GTEX_rna_region <- thisLine_GTEX[,c("gene_chr", "gene_start", "gene_end")]
      colnames(thisLine_GTEX_rna_region) <- c("chr", "start", "end")
      thisLine_GTEX_rna_region <- makeGRangesFromDataFrame(thisLine_GTEX_rna_region)
      thisSet_rna_region <- makeGRangesFromDataFrame(thisSet[,1:3])
      theIntersect <- as.data.frame(findOverlaps(thisLine_GTEX_rna_region, thisSet_rna_region))
      thisSet <- thisSet[unique(theIntersect[,2]),]

      thisLine_GTEX <- do.call("rbind", replicate(nrow(thisSet), thisLine_GTEX, simplify = FALSE))
      thisLine_GTEX$rna_region <- thisSet$rna_region
      thisLine_GTEX$tf_region <- thisSet$tf_region
      thisLine_GTEX$rna_pval <- thisSet$rna_pval
      thisLine_GTEX$numTFs <- thisSet$numTFs
      thisLine_GTEX$tf_names <- thisSet$tf_names
      thisLine_GTEX$tf_pval <- thisSet$tf_pval

      if(matVar[1,2]==thisVar[3] && patVar[1,2]==thisVar[4]) {
        thisLine_GTEX$rna_refReads <- thisSet$rna_matReads
        thisLine_GTEX$rna_altReads <- thisSet$rna_patReads
        thisLine_GTEX$tf_refReads <- thisSet$tf_matReads
        thisLine_GTEX$tf_altReads <- thisSet$tf_patReads
        d1_fullGTEX_inOurSet <- rbind(d1_fullGTEX_inOurSet, thisLine_GTEX)
      }
      if(matVar[1,2]==thisVar[4] && patVar[1,2]==thisVar[3]) {
        thisLine_GTEX$rna_refReads <- thisSet$rna_patReads
        thisLine_GTEX$rna_altReads <- thisSet$rna_matReads
        thisLine_GTEX$tf_refReads <- thisSet$tf_patReads
        thisLine_GTEX$tf_altReads <- thisSet$tf_matReads
        d1_fullGTEX_inOurSet <- rbind(d1_fullGTEX_inOurSet, thisLine_GTEX)
      }

    }
  }


}




d1_fullGTEX_inOurSet$log_rna_ratio <- log((d1_fullGTEX_inOurSet$rna_altReads+1)/(d1_fullGTEX_inOurSet$rna_refReads+1))
saveFile <- paste(outDir, "GTEX_with_rnaTFInfo_donor1_brainOnly.txt")
write.table(d1_fullGTEX_inOurSet, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
saveFile <- paste(outDir, "GTEX_with_rnaTFInfo_donor1_brainOnly.txt", sep="")
write.table(d1_fullGTEX_inOurSet, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#d1_fullGTEX_inOurSet <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

cor.test(d1_fullGTEX_inOurSet$slope, d1_fullGTEX_inOurSet$log_rna_ratio, method="pearson")
cor.test(d1_fullGTEX_inOurSet$slope, d1_fullGTEX_inOurSet$log_rna_ratio, method="spearman")


################################################################################
#Problem 1: I'm running into the issue that there are multiple slopes available
#  for each of the eQTLs, because we've got this happening in every brain region.
#Solution: For a given eQTL-gene relationship, restrict to the minimum p-value.
#
#Problem 2: I have multiple RNA regions which are associated with a given
#  eQTL. This causes problems.
#Solution: I should restrict to the single most significant RNA case, I think.
#
#Problem 3: There are still some cases where a specific RNA ratio is being
#  associated with many eQTLs of varying slopes. This is typically the result of
#  multiple eQTLs being present nearby one-another in our donor. It is unclear to
#  me what proportion of them are congruent (same sign of slope) versus competing
#  (different signs of slope).
#Solution: I should restrict to only cases where there is a single eQTL in the
#  TF-biased variant region.
#
################################################################################

d1_fullGTEX_inOurSet$uniqueNames <- paste(d1_fullGTEX_inOurSet[,1], d1_fullGTEX_inOurSet[,2], sep="_")
theList <- unique(d1_fullGTEX_inOurSet$uniqueNames)

mostSig_df_d1 <- c()
for (i in 1:length(theList)) {
  if(i%%100==0) {print(paste(i, length(theList)))}
  thisSet <- d1_fullGTEX_inOurSet[d1_fullGTEX_inOurSet$uniqueNames==theList[i],]
  thisSet <- thisSet[thisSet$pval_nominal==min(thisSet$pval_nominal),]
  thisSet <- thisSet[thisSet$rna_pval==min(thisSet$rna_pval),]
  if(length(unique(thisSet[,1]))==1) {
    mostSig_df_d1 <- rbind(mostSig_df_d1, thisSet)
  }
}



a <- mostSig_df_d1[abs(mostSig_df_d1$tss_distance)<=1000,]
cor.test(a$slope, a$log_rna_ratio)
cor.test(mostSig_df_d1$slope, mostSig_df_d1$log_rna_ratio)
cor.test(mostSig_df_d1$slope, mostSig_df_d1$log_rna_ratio, method="spearman")




################################################################################
#I would like to repeat the above, but make sure no Input significant cases
#are included.
################################################################################

pval_table_1 <- read.table(pval_table_1, header=T, sep="\t", stringsAsFactors=F)

significant_input_d1 <- pval_table_1[pval_table_1[,"INPUT"]<=0.05,]


mostSig_df_d1_noInput <- mostSig_df_d1[!mostSig_df_d1$tf_region%in%rownames(significant_input_d1),]
mostSig_df_d1_noInput <- mostSig_df_d1_noInput[!mostSig_df_d1_noInput$rna_region%in%rownames(significant_input_d1),]


a <- mostSig_df_d1_noInput[abs(mostSig_df_d1_noInput$tss_distance)<=1000,]
cor.test(a$slope, a$log_rna_ratio)
cor.test(mostSig_df_d1_noInput$slope, mostSig_df_d1_noInput$log_rna_ratio)
cor.test(mostSig_df_d1_noInput$slope, mostSig_df_d1_noInput$log_rna_ratio, method="spearman")


################################################################################
################################################################################
#Donor2
################################################################################
################################################################################


rna_tf_sig_wPhase_d2 <- read.table(rna_tf_sig_saveFile_d2, header=T, sep="\t", stringsAsFactors=F)
relevant_tf_regions_d2 <- rna_tf_sig_wPhase_d2[,4:6]
colnames(relevant_tf_regions_d2) <- c("chr", "start", "end")
relevant_tf_regions_d2 <- makeGRangesFromDataFrame(relevant_tf_regions_d2)
theIntersect <- as.data.frame(findOverlaps(GTEX_regions, relevant_tf_regions_d2))

d2_fullGTEX <- fullGTEX[unique(theIntersect[,1]),]
rna_tf_sig_wPhase_d2 <- rna_tf_sig_wPhase_d2[unique(theIntersect[,2]),]




################################################################################
#Restrict both gtex variants and our rna variants to cases where the
#RNA variation is within the appropriate gene.
################################################################################

d2_fullGTEX_rna_regions <- d2_fullGTEX[,c("gene_chr", "gene_start", "gene_end")]
colnames(d2_fullGTEX_rna_regions) <- c("chr", "start", "end")
d2_fullGTEX_rna_regions <- makeGRangesFromDataFrame(d2_fullGTEX_rna_regions)

relevant_rna_regions_d2 <- rna_tf_sig_wPhase_d2[,1:3]
colnames(relevant_rna_regions_d2) <- c("chr", "start", "end")
relevant_rna_regions_d2 <- makeGRangesFromDataFrame(relevant_rna_regions_d2)

theIntersect <- as.data.frame(findOverlaps(d2_fullGTEX_rna_regions, relevant_rna_regions_d2))

d2_fullGTEX <- d2_fullGTEX[unique(theIntersect[,1]),]
rna_tf_sig_wPhase_d2 <- rna_tf_sig_wPhase_d2[unique(theIntersect[,2]),]


################################################################################
#For ease in the next section, restrict the gtex data and
#add a tf region column to the rna_tf_sig_wPhase table
################################################################################

d2_fullGTEX <- d2_fullGTEX[,c("alt_variant_id", "gene_id", "tss_distance", "maf", "pval_nominal", "slope", "gene_name", "gene_chr", "gene_start", "gene_end", "strand", "Tissue")]

rna_tf_sig_wPhase_d2$tf_region <- paste(rna_tf_sig_wPhase_d2[,"chr_tf"], ":", rna_tf_sig_wPhase_d2[,"start_tf"], "-", rna_tf_sig_wPhase_d2[,"end_tf"], sep="")
rna_tf_sig_wPhase_d2$rna_region <- paste(rna_tf_sig_wPhase_d2[,"chr"], ":", rna_tf_sig_wPhase_d2[,"start"], "-", rna_tf_sig_wPhase_d2[,"end"], sep="")

################################################################################
#We need to restrict the table further based on actual
#variants that we observe.  To do this, I need to read in the
#names of all variant regions and parse out, for each GTEx entry:
#1) do we have a variant at the relevant location
#2) if so, do we have the appropriate reference and alternate allele?
#3) if so, correctly sort the read counts from mat/pat to ref/alt for
#   a proper comparison to GTEx predictions.
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

d2_relevant_regions <- paste(rna_tf_sig_wPhase_d2[,4], ":", rna_tf_sig_wPhase_d2[,5], "-", rna_tf_sig_wPhase_d2[,6], sep="")
d2_regions_w_vars <- d2_regions_w_vars[d2_regions_w_vars$region%in%d2_relevant_regions,]
d2_regions_w_vars_gr <- makeGRangesFromDataFrame(d2_regions_w_vars)



d2_fullGTEX_inOurSet <- c()

for (i in 1:nrow(d2_fullGTEX)) {
  if(i%%100==0) {print(paste(i, nrow(d2_fullGTEX), nrow(d2_fullGTEX_inOurSet)))}
  thisVar <- strsplit(d2_fullGTEX[i,"alt_variant_id"], split="_")[[1]]
  thisVar_gr <- data.frame(chr=thisVar[1], start=as.numeric(thisVar[2]), end=as.numeric(thisVar[2]))
  thisVar_gr <- makeGRangesFromDataFrame(thisVar_gr)

  theIntersect <- as.data.frame(findOverlaps(d2_regions_w_vars_gr, thisVar_gr))
  if(nrow(theIntersect)>0) {
    tf_region <- d2_regions_w_vars[theIntersect[1,1],"region"]

    matVar <- d2_regions_w_vars[theIntersect[1,1],"matVar"]
    patVar <- d2_regions_w_vars[theIntersect[1,1],"patVar"]

    matVar <- strsplit(matVar, split=";")[[1]]
    patVar <- strsplit(patVar, split=";")[[1]]

    matVar <- matrix(nrow=length(matVar), ncol=3, data=unlist(strsplit(matVar, split=",")), byrow=T)
    patVar <- matrix(nrow=length(patVar), ncol=3, data=unlist(strsplit(patVar, split=",")), byrow=T)

    matVar <- as.data.frame(matVar)
    colnames(matVar) <- c("loc", "nuc", "type")
    patVar <- as.data.frame(patVar)
    colnames(patVar) <- c("loc", "nuc", "type")

    matVar <- matVar[matVar[,3]=="HET",]
    patVar <- patVar[patVar[,3]=="HET",]

    matVar <- matVar[matVar[,"loc"]==thisVar[2],]
    patVar <- patVar[patVar[,"loc"]==thisVar[2],]


    if(nrow(matVar)==1 && nrow(patVar)==1) {
      thisLine_GTEX <- d2_fullGTEX[i,]
      thisSet <- rna_tf_sig_wPhase_d2[rna_tf_sig_wPhase_d2$tf_region==tf_region,]

      thisLine_GTEX_rna_region <- thisLine_GTEX[,c("gene_chr", "gene_start", "gene_end")]
      colnames(thisLine_GTEX_rna_region) <- c("chr", "start", "end")
      thisLine_GTEX_rna_region <- makeGRangesFromDataFrame(thisLine_GTEX_rna_region)
      thisSet_rna_region <- makeGRangesFromDataFrame(thisSet[,1:3])
      theIntersect <- as.data.frame(findOverlaps(thisLine_GTEX_rna_region, thisSet_rna_region))
      thisSet <- thisSet[unique(theIntersect[,2]),]

      thisLine_GTEX <- do.call("rbind", replicate(nrow(thisSet), thisLine_GTEX, simplify = FALSE))
      thisLine_GTEX$rna_region <- thisSet$rna_region
      thisLine_GTEX$tf_region <- thisSet$tf_region
      thisLine_GTEX$rna_pval <- thisSet$rna_pval
      thisLine_GTEX$numTFs <- thisSet$numTFs
      thisLine_GTEX$tf_names <- thisSet$tf_names
      thisLine_GTEX$tf_pval <- thisSet$tf_pval

      if(matVar[1,2]==thisVar[3] && patVar[1,2]==thisVar[4]) {
        thisLine_GTEX$rna_refReads <- thisSet$rna_matReads
        thisLine_GTEX$rna_altReads <- thisSet$rna_patReads
        thisLine_GTEX$tf_refReads <- thisSet$tf_matReads
        thisLine_GTEX$tf_altReads <- thisSet$tf_patReads
        d2_fullGTEX_inOurSet <- rbind(d2_fullGTEX_inOurSet, thisLine_GTEX)
      }
      if(matVar[1,2]==thisVar[4] && patVar[1,2]==thisVar[3]) {
        thisLine_GTEX$rna_refReads <- thisSet$rna_patReads
        thisLine_GTEX$rna_altReads <- thisSet$rna_matReads
        thisLine_GTEX$tf_refReads <- thisSet$tf_patReads
        thisLine_GTEX$tf_altReads <- thisSet$tf_matReads
        d2_fullGTEX_inOurSet <- rbind(d2_fullGTEX_inOurSet, thisLine_GTEX)
      }

    }
  }


}




d2_fullGTEX_inOurSet$log_rna_ratio <- log((d2_fullGTEX_inOurSet$rna_altReads+1)/(d2_fullGTEX_inOurSet$rna_refReads+1))
saveFile <- paste(outDir, "GTEX_with_rnaTFInfo_donor2_brainOnly.txt")
write.table(d2_fullGTEX_inOurSet, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
saveFile <- paste(outDir, "GTEX_with_rnaTFInfo_donor2_brainOnly.txt", sep="")
write.table(d2_fullGTEX_inOurSet, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#d2_fullGTEX_inOurSet <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)


cor.test(d2_fullGTEX_inOurSet$slope, d2_fullGTEX_inOurSet$log_rna_ratio, method="pearson")
cor.test(d2_fullGTEX_inOurSet$slope, d2_fullGTEX_inOurSet$log_rna_ratio, method="spearman")

################################################################################
#Problem 1: I'm running into the issue that there are multiple slopes available
#  for each of the eQTLs, because we've got this happening in every brain region.
#Solution: For a given eQTL-gene relationship, restrict to the minimum p-value.
#
#Problem 2: I have multiple RNA regions which are associated with a given
#  eQTL. This causes problems.
#Solution: I should restrict to the single most significant RNA case, I think.
#
#Problem 3: There are still some cases where a specific RNA ratio is being
#  associated with many eQTLs of varying slopes. This is typically the result of
#  multiple eQTLs being present nearby one-another in our donor. It is unclear to
#  me what proportion of them are congruent (same sign of slope) versus competing
#  (different signs of slope).
#Solution: I should restrict to only cases where there is a single eQTL in the
#  TF-biased variant region.
#
################################################################################

d2_fullGTEX_inOurSet$uniqueNames <- paste(d2_fullGTEX_inOurSet[,1], d2_fullGTEX_inOurSet[,2], sep="_")
theList <- unique(d2_fullGTEX_inOurSet$uniqueNames)

mostSig_df_d2 <- c()
for (i in 1:length(theList)) {
  if(i%%100==0) {print(paste(i, length(theList)))}
  thisSet <- d2_fullGTEX_inOurSet[d2_fullGTEX_inOurSet$uniqueNames==theList[i],]
  thisSet <- thisSet[thisSet$pval_nominal==min(thisSet$pval_nominal),]
  thisSet <- thisSet[thisSet$rna_pval==min(thisSet$rna_pval),]
  if(length(unique(thisSet[,1]))==1) {
    mostSig_df_d2 <- rbind(mostSig_df_d2, thisSet)
  }
}


a <- mostSig_df_d2[abs(mostSig_df_d2$tss_distance)<=1000,]
cor.test(a$slope, a$log_rna_ratio)
cor.test(mostSig_df_d2$slope, mostSig_df_d2$log_rna_ratio)
cor.test(mostSig_df_d2$slope, mostSig_df_d2$log_rna_ratio, method="spearman")



################################################################################
#I would like to repeat the above, but make sure no Input significant cases
#are included.
################################################################################

pval_table_2 <- read.table(pval_table_2, header=T, sep="\t", stringsAsFactors=F)

significant_input_d2 <- pval_table_2[pval_table_2[,"INPUT"]<=0.05,]


mostSig_df_d2_noInput <- mostSig_df_d2[!mostSig_df_d2$tf_region%in%rownames(significant_input_d2),]
mostSig_df_d2_noInput <- mostSig_df_d2_noInput[!mostSig_df_d2_noInput$rna_region%in%rownames(significant_input_d2),]


a <- mostSig_df_d2_noInput[abs(mostSig_df_d2_noInput$tss_distance)<=1000,]
cor.test(a$slope, a$log_rna_ratio)
cor.test(mostSig_df_d2_noInput$slope, mostSig_df_d2_noInput$log_rna_ratio)
cor.test(mostSig_df_d2_noInput$slope, mostSig_df_d2_noInput$log_rna_ratio, method="spearman")





################################################################################
################################################################################
#Shared and combination.
################################################################################
################################################################################


################################################################################
#Greg has requested a rework of the above figure to be a binned violin plot.
#Make sure no Input significant cases are included.
################################################################################

################################################################################
#I would like to repeat the above, but make sure no Input significant cases
#are included.
################################################################################




graphDF <- rbind(mostSig_df_d1_noInput[,c("slope", "log_rna_ratio")], mostSig_df_d2_noInput[,c("slope", "log_rna_ratio")])

#theBins <- c(-2, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 2)
theBins <- c(0.05, 0.25, 0.5, 0.75, 1)
graphDF$slopeBin <- rep(0, nrow(graphDF))

for (i in 1:length(theBins)) {
  graphDF[graphDF$slope< -theBins[i],"slopeBin"] <- -theBins[i]
  graphDF[graphDF$slope> theBins[i],"slopeBin"] <- theBins[i]
}

for (i in 1:length(theBins)) {
  graphDF[graphDF$slope>theBins[i],"slopeBin"] <- theBins[i]
}


table(graphDF$slopeBin)

graphDF$slopeBin <- factor(graphDF$slopeBin, levels=c(-1, -0.75, -0.5, -0.25, -0.05, 0.05, 0.25, 0.5, 0.75, 1))

#saveFile <- paste(outDir, "GTEX_v_rnaRatio_violinBins_combined_brainOnly_noInput.pdf", sep="")
#p <- ggplot(graphDF, aes(y=log_rna_ratio, x=slopeBin)) + geom_violin(alpha=1) + theme_classic(base_size = 20) +
#  ylim(-7.5,7.5) + xlab("GTEx Slope") + ylab ("log(RNA Bias)") + theme(legend.position = "none")
#ggsave(saveFile)


saveFile <- paste(outDir, "Figure_4A.pdf", sep="")
p <- ggplot(graphDF, aes(y=log_rna_ratio, x=slopeBin)) + geom_boxplot() + theme_classic(base_size = 20) +
  ylim(-7.5,7.5) + xlab("GTEx Slope") + ylab ("log(RNA Bias)") + theme(legend.position = "none")
ggsave(saveFile)




#
