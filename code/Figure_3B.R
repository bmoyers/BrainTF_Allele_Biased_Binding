#!/usr/bin/env R
#Figure_3B.R

################################################################################
#This script is used to produce Figures 3B.
#
#Run under R 4.1.0, but can be used under other versions so long as the packages
#     GenomicRanges, matrixStats, and ggplot are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# vcf_dir: Directory containing the vcf files from ensembl at the following URL:
#     https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/ ,
#     in August of 2022.
# donor_variant_table_d1_saveFile: the VCF file for donor 1.
# pval_table_1: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# depthTable_1_mat: A table compiling all TF read  depths for donor 1, haplotype 1,
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depthTable_1_pat: As depthTable_1_mat, but for haplotype 2.
# vep_d1_saveFile: a vep-annotated vcf file for donor 1, produced by vep using the
#     vep108.ini file.
# allSeqs_d1_saveFile: path to the fasta file for all variant sequences in donor 1,
#     produced by the script Building_fasta_variant_sequences.sh
# donor_variant_table_d2_saveFile: the VCF file for donor 2.
# pval_table_2: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 2, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# depthTable_2_mat: A table compiling all TF read  depths for donor 2, haplotype 1,
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depthTable_2_pat: As depthTable_1_mat, but for haplotype 2.
# vep_d2_saveFile: a vep-annotated vcf file for donor 2, produced by vep using the
#     vep108.ini file.
# allSeqs_d2_saveFile: path to the fasta file for all variant sequences in donor 2,
#     produced by the script Building_fasta_variant_sequences.sh
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



convert_regions_to_gr <- function(peakFile, totalCol=3, theNames=c("chr", "start", "end")) {
  peakFile <- as.data.frame(peakFile[,1:totalCol])
  peakFile[,1] <- as.character(peakFile[,1])
  peakFile[,2] <- as.numeric(as.character(peakFile[,2]))
  peakFile[,3] <- as.numeric(as.character(peakFile[,3]))
  colnames(peakFile) <- theNames
  peakFile_gr <- makeGRangesFromDataFrame(peakFile, keep.extra.columns=T, ignore.strand=T)
  return(peakFile_gr)
}



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################



args <- commandArgs(trailingOnly=T)
outDir <- args[1]
vcf_dir <- args[2]
donor_variant_table_d1 <- args[3]
pval_table_1 <- args[4]
depthTable_1_mat <- args[5]
depthTable_1_pat <- args[6]
vep_d1_saveFile <- args[7]
allSeqs_d1_saveFile <- args[8]
donor_variant_table_d2 <- args[9]
pval_table_2 <- args[10]
depthTable_2_mat <- args[11]
depthTable_2_pat <- args[12]
vep_d2_saveFile <- args[13]
allSeqs_d2_saveFile <- args[14]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#vcf_dir <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/ensembl/"
#donor_variant_table_d1_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants.vcf.table"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#depthTable_1_mat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depthTable_1_pat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#vep_d1_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0001_annotated.tsv.gz"
#allSeqs_d1_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#donor_variant_table_d2_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.table"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#depthTable_2_mat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depthTable_2_pat <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#vep_d2_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0002_annotated.tsv.gz"
#allSeqs_d2_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"


################################################################################
#For each donor, remove histone marks and pol.
#identify the minimum significance level of each variant.
#identify the DAF of the variant/variants within that region.
#Create a table and barplot of significance as a function of DAF.
#Make this a frequency plot among only the significant cases for simplicity.
################################################################################



################################################################################
#Get VCF information for variants near our donor variants.
################################################################################


################################################################################
#Note that this first section need not be repeated-- I was grabbing the
#newest ancestral allele information and preparing it so I can read things
#in efficiently.
################################################################################

donor_variant_table_d1 <- read.table(donor_variant_table_d1_saveFile, header=F, sep="\t", stringsAsFactors=F)
donor_variant_table_d1$name <- paste(donor_variant_table_d1[,1], donor_variant_table_d1[,2], sep="_")



donor_variant_table_d2 <- read.table(donor_variant_table_d2_saveFile, header=F, sep="\t", stringsAsFactors=F)
donor_variant_table_d2$name <- paste(donor_variant_table_d2[,1], donor_variant_table_d2[,2], sep="_")

donor_variant_table_both <- rbind(donor_variant_table_d1, donor_variant_table_d2)

uniqueChroms <- unique(donor_variant_table_both[,1])
uniqueChroms <- uniqueChroms[!uniqueChroms%in%c("chrY", "chrX")]


vcf_in_our_donors_with_ancestral_allele <- c()
for (i in 1:length(uniqueChroms)) {
  print(paste(i, length(uniqueChroms), uniqueChroms[i], nrow(vcf_in_our_donors_with_ancestral_allele)))
  thisVCF <- paste(vcf_dir, "homo_sapiens-", uniqueChroms[i], ".vcf.gz", sep="")
  thisVCF <- read.table(thisVCF, comment.char="#", header=F, sep="\t", stringsAsFactors=F)
  colnames(thisVCF) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  thisVCF[,1] <- paste("chr", thisVCF[,1], sep="")
  thisVCF$name <- paste(thisVCF[,1], thisVCF[,2], sep="_")
  thisVCF <- thisVCF[thisVCF$name%in%donor_variant_table_both$name,]
  thisVCF  <- thisVCF[grep("AA=", thisVCF[,"INFO"]),]
  vcf_in_our_donors_with_ancestral_allele <- rbind(vcf_in_our_donors_with_ancestral_allele, thisVCF)
}


saveFile <- paste(vcf_dir, "vcf_all_variants_in_donors_1_and_2_with_AncestralAllele_info.txt", sep="")
write.table(vcf_in_our_donors_with_ancestral_allele, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#vcf_in_our_donors_with_ancestral_allele <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)



################################################################################
#Read in the specific variants in our donors.
################################################################################

donor_variant_table_d1 <- read.table(donor_variant_table_d1_saveFile, header=F, sep="\t", stringsAsFactors=F)
donor_variant_table_d1$name <- paste(donor_variant_table_d1[,1], donor_variant_table_d1[,2], sep="_")

donor_variant_table_d2 <- read.table(donor_variant_table_d2_saveFile, header=F, sep="\t", stringsAsFactors=F)
donor_variant_table_d2$name <- paste(donor_variant_table_d2[,1], donor_variant_table_d2[,2], sep="_")

donor_variant_table_both <- rbind(donor_variant_table_d1, donor_variant_table_d2)


vcf_in_our_donors_with_ancestral_allele_df <- vcf_in_our_donors_with_ancestral_allele[,c(1,2,2)]
colnames(vcf_in_our_donors_with_ancestral_allele_df) <- c("chr", "start", "end")
vcf_in_our_donors_with_ancestral_allele_gr <- makeGRangesFromDataFrame(vcf_in_our_donors_with_ancestral_allele_df)


AncAll <- rep(NA, nrow(vcf_in_our_donors_with_ancestral_allele))
for (i in 1:nrow(vcf_in_our_donors_with_ancestral_allele)) {
  if(i%%1000==0) { print(paste(i, nrow(vcf_in_our_donors_with_ancestral_allele)))}
  thisInfo <- strsplit(vcf_in_our_donors_with_ancestral_allele[i,"INFO"], split=";")[[1]]
  thisAncAll <- thisInfo[grep("AA=", thisInfo)]
  thisAncAll <- gsub("AA=", "", thisAncAll)
  AncAll[i] <- thisAncAll
}

vcf_in_our_donors_with_ancestral_allele$AncAll <- AncAll




################################################################################
################################################################################
#Construct a table with all relevant info for Donor 1.
################################################################################
################################################################################

pval_table_1 <- read.table(pval_table_1, header=T, sep="\t", stringsAsFactors=F)
depthTable_1_mat <- read.table(depthTable_1_mat, header=T, sep="\t", stringsAsFactors=F)
depthTable_1_pat <- read.table(depthTable_1_pat, header=T, sep="\t", stringsAsFactors=F)

#Identify cases where INPUT was significant

significant_input_d1 <- pval_table_1[pval_table_1[,"INPUT"]<=0.05,]
significant_input_d1_stricter <- pval_table_1[pval_table_1[,"INPUT"]<=0.001,]


toRemove <- grep("chrX", rownames(pval_table_1))
if(length(toRemove)>0) {
  pval_table_1 <- pval_table_1[-toRemove,]
  depthTable_1_mat <- depthTable_1_mat[-toRemove,]
  depthTable_1_pat <- depthTable_1_pat[-toRemove,]
}


forbiddenExprs_2 <- c(colnames(pval_table_1)[grep("H[3|4]K", colnames(pval_table_1))], colnames(pval_table_1)[grep("INPUT", colnames(pval_table_1))], colnames(pval_table_1)[grep("POL2", colnames(pval_table_1))])
pval_table_1 <- pval_table_1[,!colnames(pval_table_1)%in%forbiddenExprs_2]


depthTable_1_mat <- depthTable_1_mat[,colnames(depthTable_1_mat)%in%colnames(pval_table_1)]
depthTable_1_pat <- depthTable_1_pat[,colnames(depthTable_1_pat)%in%colnames(pval_table_1)]


depthTable_1_mat <- depthTable_1_mat[rowMins(as.matrix(pval_table_1))<=0.05,]
depthTable_1_pat <- depthTable_1_pat[rowMins(as.matrix(pval_table_1))<=0.05,]
pval_table_1 <- pval_table_1[rowMins(as.matrix(pval_table_1))<=0.05,]

depthTable_1 <- depthTable_1_mat + depthTable_1_pat



all_significance_effects <- c()

for (i in 1:nrow(pval_table_1)) {
  if(i%%100==0) {print(paste(i, nrow(pval_table_1)))}
  minPval_options <- which(as.numeric(pval_table_1[i,])==min(as.numeric(pval_table_1[i,])))

  this_Table <- cbind(as.numeric(pval_table_1[i,minPval_options]), as.numeric(depthTable_1_mat[i,minPval_options]), as.numeric(depthTable_1_pat[i,minPval_options]), as.numeric(depthTable_1[i,minPval_options]))
  this_Table <- as.data.frame(this_Table)
  colnames(this_Table) <- c("these_pvals", "matDepth", "patDepth", "totalDepth")
  this_Table$Regions <- rep(rownames(pval_table_1)[i], nrow(this_Table))
  for (j in 1:(ncol(this_Table)-1)) { this_Table[,j] <- as.numeric(as.character(this_Table[,j])) }

  this_Table$effectSize <- rowMaxs(as.matrix(this_Table[,2:3]))/this_Table$totalDepth
  this_Table$Donor <- rep("Donor1", nrow(this_Table))
  this_Table <- unique(this_Table)
  if(nrow(this_Table)>1) {
    print(i)
    #break()
  }
  all_significance_effects <- rbind(all_significance_effects, this_Table)

}



all_significance_effects_regions <- matrix(nrow=nrow(all_significance_effects), ncol=3, data=unlist(strsplit(all_significance_effects$Regions, split=":|-")), byrow=T)
all_significance_effects_regions_gr <- convert_regions_to_gr(all_significance_effects_regions)

theIntersect <- as.data.frame(findOverlaps(all_significance_effects_regions_gr, vcf_in_our_donors_with_ancestral_allele_gr))

all_significance_effects_new <- cbind(all_significance_effects[theIntersect[,1],], vcf_in_our_donors_with_ancestral_allele[theIntersect[,2],c("CHROM", "POS", "REF", "ALT", "AncAll")])

sine_matches <- all_significance_effects[-theIntersect[,1],]

sine_matches$CHROM <- rep(NA, nrow(sine_matches))
sine_matches$POS <- rep(NA, nrow(sine_matches))
sine_matches$REF <- rep(NA, nrow(sine_matches))
sine_matches$ALT <- rep(NA, nrow(sine_matches))
sine_matches$AncAll <- rep(NA, nrow(sine_matches))


all_significance_effects_d1 <- rbind(all_significance_effects_new, sine_matches)


rm(donor_variant_table_d1, AncAll, all_significance_effects_new, all_significance_effects_regions_gr, depthTable_1_pat, sine_matches, all_significance_effects, all_significance_effects_regions, depthTable_1_mat, pval_table_1)

################################################################################
#Determine which allele was assigned Mat and which was assigned Pat.
################################################################################

allSeqs_d1 <- readLines(allSeqs_d1_saveFile)
allSeqs_d1 <- allSeqs_d1[grep(">", allSeqs_d1)]
allSeqs_d1 <- gsub(">", "", allSeqs_d1)
allSeqs_d1_mat <- allSeqs_d1[grep(":M:", allSeqs_d1, fixed=T)]
allSeqs_d1_pat <- allSeqs_d1[grep(":P:", allSeqs_d1, fixed=T)]

allSeqs_d1_regions <- matrix(nrow=length(allSeqs_d1_mat), ncol=2, data=unlist(strsplit(allSeqs_d1_mat, split=":M:")), byrow=T)
allSeqs_d1_regions <- matrix(nrow=nrow(allSeqs_d1_regions), ncol=3, data=unlist(strsplit(allSeqs_d1_regions[,1], split=":|-")), byrow=T)

allSeqs_d1_mat <- matrix(nrow=length(allSeqs_d1_mat), ncol=2, data=unlist(strsplit(allSeqs_d1_mat, split=":M:")), byrow=T)
allSeqs_d1_pat <- matrix(nrow=length(allSeqs_d1_pat), ncol=2, data=unlist(strsplit(allSeqs_d1_pat, split=":P:")), byrow=T)


allSeqs_d1_allInfo <- cbind(allSeqs_d1_regions, allSeqs_d1_mat, allSeqs_d1_pat[,2])
allSeqs_d1_allInfo <- as.data.frame(allSeqs_d1_allInfo)
colnames(allSeqs_d1_allInfo) <- c("chr", "start", "end", "name", "Mat", "Pat")

rm(allSeqs_d1, allSeqs_d1_mat, allSeqs_d1_pat, allSeqs_d1_regions)



matAllele <- rep(NA, nrow(all_significance_effects_d1))
patAllele <- rep(NA, nrow(all_significance_effects_d1))

theMatch <- match(all_significance_effects_d1$Regions, allSeqs_d1_allInfo$name)

theMats <- as.character(allSeqs_d1_allInfo[theMatch,"Mat"])
thePats <- as.character(allSeqs_d1_allInfo[theMatch,"Pat"])

problemChildren <- c()
for (i in 1:length(theMats)) {
  if(i%%1000==0) {print(paste(i, length(theMats)))}
  thisMatAllele <- strsplit(theMats[i], split=";")[[1]]
  thisPatAllele <- strsplit(thePats[i], split=";")[[1]]

  thisMatAllele <- matrix(nrow=length(thisMatAllele), ncol=3, data=unlist(strsplit(thisMatAllele, split=",")), byrow=T)
  thisPatAllele <- matrix(nrow=length(thisPatAllele), ncol=3, data=unlist(strsplit(thisPatAllele, split=",")), byrow=T)

  thisLoc <- all_significance_effects_d1[i,"POS"]
  if(!thisLoc%in%thisMatAllele[,1]) { problemChildren <- c(problemChildren, i)}
  if(thisLoc%in%thisMatAllele[,1]) {
    matAllele[i] <- thisMatAllele[thisMatAllele[,1]==all_significance_effects_d1[i,"POS"],2]
    patAllele[i] <- thisPatAllele[thisPatAllele[,1]==all_significance_effects_d1[i,"POS"],2]
  }
}

all_significance_effects_d1$MatAllele <- matAllele
all_significance_effects_d1$PatAllele <- patAllele


################################################################################
#Now we add the allele frequency, where data was available.
#Note that this needs to be corrected-- while we identify the Allele Frequency
#via gnomAD, we have to then correct for which of the two alleles is
#ancestral versus derived.
################################################################################

vep_d1 <- readLines(vep_d1_saveFile)
vep_d1 <- vep_d1[-grep("^##", vep_d1)]
vep_d1_header <- vep_d1[grep("^#", vep_d1)]
vep_d1_header <- gsub("^#", "", vep_d1_header)
vep_d1_header <- strsplit(vep_d1_header, split="\t")[[1]]

vep_d1 <- read.delim(vep_d1_saveFile, header=F, sep="\t", stringsAsFactors=F, comment.char="#")
colnames(vep_d1) <- vep_d1_header

vep_d1 <- vep_d1[,c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
      "Existing_variation", "DISTANCE", "STRAND", "VARIANT_CLASS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE",
      "CANONICAL", "gnomADg_AF", "CLIN_SIG", "PUBMED", "CADD_PHRED", "CADD_RAW", "gnomad3_AF", "GerpRS")]


all_significance_effects_d1$name1 <- paste(all_significance_effects_d1$CHROM, "_", all_significance_effects_d1$POS, "_", all_significance_effects_d1$MatAllele, "/", all_significance_effects_d1$PatAllele, sep="")
all_significance_effects_d1$name2 <- paste(all_significance_effects_d1$CHROM, "_", all_significance_effects_d1$POS, "_", all_significance_effects_d1$PatAllele, "/", all_significance_effects_d1$MatAllele, sep="")


match1 <- match(all_significance_effects_d1$name1, vep_d1$Uploaded_variation)
match2 <- match(all_significance_effects_d1$name2, vep_d1$Uploaded_variation)



Alt_AF1 <- rep(NA, nrow(all_significance_effects_d1))
Alt_AF1 <- vep_d1[match1, "gnomad3_AF"]

Alt_AF2 <- rep(NA, nrow(all_significance_effects_d1))
Alt_AF2 <- vep_d1[match2, "gnomad3_AF"]

Alt_AF <- cbind(Alt_AF1, Alt_AF2)
Alt_AF <- as.data.frame(Alt_AF)
Alt_AF[,1] <- as.numeric(as.character(Alt_AF[,1]))
Alt_AF[,2] <- as.numeric(as.character(Alt_AF[,2]))

final_Alt_AF <- rep(NA, nrow(Alt_AF))
for (i in 1:nrow(Alt_AF)) {
  if(!is.na(Alt_AF[i,1]) && is.na(Alt_AF[i,2])) { final_Alt_AF[i] <- Alt_AF[i,1]}
  if(is.na(Alt_AF[i,1]) && !is.na(Alt_AF[i,2])) { final_Alt_AF[i] <- Alt_AF[i,2]}
  if(!is.na(Alt_AF[i,1]) && !is.na(Alt_AF[i,2])) {
    print(i)
    break()
  }

}


all_significance_effects_d1$Alt_AF <- final_Alt_AF



name1 <- rep(NA, nrow(all_significance_effects_d1))
name2 <- rep(NA, nrow(all_significance_effects_d1))
for (i in 1:nrow(all_significance_effects_d1)) {
  if(i%%1000==0) {print(paste(i, nrow(all_significance_effects_d1)))}
  twoOptions <- all_significance_effects_d1[i,c("MatAllele", "PatAllele")]
  theAncestralAllele <- all_significance_effects_d1[i,"AncAll"]
  theDerivedAllele <- twoOptions[twoOptions!=theAncestralAllele]
  if(length(theDerivedAllele)==1) {
    this_name1 <- paste(all_significance_effects_d1[i,"CHROM"], "_", all_significance_effects_d1[i,"POS"], "_", theAncestralAllele, "/", theDerivedAllele, sep="")
    name1[i] <- this_name1
    this_name2 <- paste(all_significance_effects_d1[i,"CHROM"], "_", all_significance_effects_d1[i,"POS"], "_", theDerivedAllele, "/", theAncestralAllele, sep="")
    name2[i] <- this_name2
  }

}



all_significance_effects_d1$name1 <- name1
all_significance_effects_d1$name2 <- name2


match1 <- match(all_significance_effects_d1$name1, vep_d1$Uploaded_variation)
match2 <- match(all_significance_effects_d1$name2, vep_d1$Uploaded_variation)


DAF1 <- rep(NA, nrow(all_significance_effects_d1))
DAF1 <- vep_d1[match1, "gnomad3_AF"]

DAF2 <- rep(NA, nrow(all_significance_effects_d1))
DAF2 <- vep_d1[match2, "gnomad3_AF"]

DAF <- cbind(DAF1, DAF2)
DAF <- as.data.frame(DAF)
DAF[,1] <- as.numeric(as.character(DAF[,1]))
DAF[,2] <- as.numeric(as.character(DAF[,2]))

#When the derived is the reference, we need to subtract from 1.
DAF[,2] <- 1-DAF[,2]


finalDAF <- rep(NA, nrow(DAF))
for (i in 1:nrow(DAF)) {
  if(!is.na(DAF[i,1]) && is.na(DAF[i,2])) { finalDAF[i] <- DAF[i,1]}
  if(is.na(DAF[i,1]) && !is.na(DAF[i,2])) { finalDAF[i] <- DAF[i,2]}
  if(!is.na(DAF[i,1]) && !is.na(DAF[i,2])) {
    print(i)
    break()
  }

}


all_significance_effects_d1$DAF <- finalDAF

rm(DAF, DAF1, DAF2, finalDAF, matAllele, patAllele, match1, match2, theMatch, theMats, thePats, vep_d1, vep_d1_header)


all_significance_effects_d1 <- all_significance_effects_d1[!is.na(all_significance_effects_d1$MatAllele),]
all_significance_effects_d1 <- all_significance_effects_d1[!is.na(all_significance_effects_d1$PatAllele),]



################################################################################
################################################################################
#Construct a table with all relevant info for Donor 2.
################################################################################
################################################################################

pval_table_2 <- read.table(pval_table_2, header=T, sep="\t", stringsAsFactors=F)
depthTable_2_mat <- read.table(depthTable_2_mat, header=T, sep="\t", stringsAsFactors=F)
depthTable_2_pat <- read.table(depthTable_2_pat, header=T, sep="\t", stringsAsFactors=F)


#Identify cases where INPUT was significant

significant_input_d2 <- pval_table_2[pval_table_2[,"INPUT"]<=0.05,]
significant_input_d2_stricter <- pval_table_2[pval_table_2[,"INPUT"]<=0.001,]


toRemove <- grep("chrX", rownames(pval_table_2))
if(length(toRemove)>0) {
  pval_table_2 <- pval_table_2[-toRemove,]
  depthTable_2_mat <- depthTable_2_mat[-toRemove,]
  depthTable_2_pat <- depthTable_2_pat[-toRemove,]
}



forbiddenExprs_2 <- c(colnames(pval_table_2)[grep("H3K", colnames(pval_table_2))], colnames(pval_table_2)[grep("INPUT", colnames(pval_table_2))], colnames(pval_table_2)[grep("POL2", colnames(pval_table_2))])
pval_table_2 <- pval_table_2[,!colnames(pval_table_2)%in%forbiddenExprs_2]


depthTable_2_mat <- depthTable_2_mat[,colnames(depthTable_2_mat)%in%colnames(pval_table_2)]
depthTable_2_pat <- depthTable_2_pat[,colnames(depthTable_2_pat)%in%colnames(pval_table_2)]



depthTable_2_mat <- depthTable_2_mat[rowMins(as.matrix(pval_table_2))<=0.05,]
depthTable_2_pat <- depthTable_2_pat[rowMins(as.matrix(pval_table_2))<=0.05,]
pval_table_2 <- pval_table_2[rowMins(as.matrix(pval_table_2))<=0.05,]

depthTable_2 <- depthTable_2_mat + depthTable_2_pat



all_significance_effects <- c()

for (i in 1:nrow(pval_table_2)) {
  if(i%%100==0) {print(paste(i, nrow(pval_table_2)))}
  minPval_options <- which(as.numeric(pval_table_2[i,])==min(as.numeric(pval_table_2[i,])))

  this_Table <- cbind(as.numeric(pval_table_2[i,minPval_options]), as.numeric(depthTable_2_mat[i,minPval_options]), as.numeric(depthTable_2_pat[i,minPval_options]), as.numeric(depthTable_2[i,minPval_options]))
  this_Table <- as.data.frame(this_Table)
  colnames(this_Table) <- c("these_pvals", "matDepth", "patDepth", "totalDepth")
  this_Table$Regions <- rep(rownames(pval_table_2)[i], nrow(this_Table))
  for (j in 1:(ncol(this_Table)-1)) { this_Table[,j] <- as.numeric(as.character(this_Table[,j])) }

  this_Table$effectSize <- rowMaxs(as.matrix(this_Table[,2:3]))/this_Table$totalDepth
  this_Table$Donor <- rep("Donor2", nrow(this_Table))
  this_Table <- unique(this_Table)
  if(nrow(this_Table)>1) {
    print(i)
    #break()
  }
  all_significance_effects <- rbind(all_significance_effects, this_Table)

}




all_significance_effects_regions <- matrix(nrow=nrow(all_significance_effects), ncol=3, data=unlist(strsplit(all_significance_effects$Regions, split=":|-")), byrow=T)
all_significance_effects_regions_gr <- convert_regions_to_gr(all_significance_effects_regions)

theIntersect <- as.data.frame(findOverlaps(all_significance_effects_regions_gr, vcf_in_our_donors_with_ancestral_allele_gr))

all_significance_effects_new <- cbind(all_significance_effects[theIntersect[,1],], vcf_in_our_donors_with_ancestral_allele[theIntersect[,2],c("CHROM", "POS", "REF", "ALT", "AncAll")])

sine_matches <- all_significance_effects[-theIntersect[,1],]

sine_matches$CHROM <- rep(NA, nrow(sine_matches))
sine_matches$POS <- rep(NA, nrow(sine_matches))
sine_matches$REF <- rep(NA, nrow(sine_matches))
sine_matches$ALT <- rep(NA, nrow(sine_matches))
sine_matches$AncAll <- rep(NA, nrow(sine_matches))


all_significance_effects_d2 <- rbind(all_significance_effects_new, sine_matches)


rm(all_significance_effects_regions_gr, depthTable_2_pat, sine_matches, all_significance_effects, all_significance_effects_regions, depthTable_2_mat, pval_table_2)



################################################################################
#Determine which allele was assigned Mat and which was assigned Pat.
################################################################################


allSeqs_d2 <- readLines(allSeqs_d2_saveFile)
allSeqs_d2 <- allSeqs_d2[grep(">", allSeqs_d2)]
allSeqs_d2 <- gsub(">", "", allSeqs_d2)
allSeqs_d2_mat <- allSeqs_d2[grep(":M:", allSeqs_d2, fixed=T)]
allSeqs_d2_pat <- allSeqs_d2[grep(":P:", allSeqs_d2, fixed=T)]

allSeqs_d2_regions <- matrix(nrow=length(allSeqs_d2_mat), ncol=2, data=unlist(strsplit(allSeqs_d2_mat, split=":M:")), byrow=T)
allSeqs_d2_regions <- matrix(nrow=nrow(allSeqs_d2_regions), ncol=3, data=unlist(strsplit(allSeqs_d2_regions[,1], split=":|-")), byrow=T)

allSeqs_d2_mat <- matrix(nrow=length(allSeqs_d2_mat), ncol=2, data=unlist(strsplit(allSeqs_d2_mat, split=":M:")), byrow=T)
allSeqs_d2_pat <- matrix(nrow=length(allSeqs_d2_pat), ncol=2, data=unlist(strsplit(allSeqs_d2_pat, split=":P:")), byrow=T)


allSeqs_d2_allInfo <- cbind(allSeqs_d2_regions, allSeqs_d2_mat, allSeqs_d2_pat[,2])
allSeqs_d2_allInfo <- as.data.frame(allSeqs_d2_allInfo)
colnames(allSeqs_d2_allInfo) <- c("chr", "start", "end", "name", "Mat", "Pat")

rm(allSeqs_d2, allSeqs_d2_mat, allSeqs_d2_pat, allSeqs_d2_regions)



matAllele <- rep(NA, nrow(all_significance_effects_d2))
patAllele <- rep(NA, nrow(all_significance_effects_d2))

theMatch <- match(all_significance_effects_d2$Regions, allSeqs_d2_allInfo$name)

theMats <- as.character(allSeqs_d2_allInfo[theMatch,"Mat"])
thePats <- as.character(allSeqs_d2_allInfo[theMatch,"Pat"])

problemChildren <- c()
for (i in 1:length(theMats)) {
  if(i%%1000==0) {print(paste(i, length(theMats)))}
  thisMatAllele <- strsplit(theMats[i], split=";")[[1]]
  thisPatAllele <- strsplit(thePats[i], split=";")[[1]]

  thisMatAllele <- matrix(nrow=length(thisMatAllele), ncol=3, data=unlist(strsplit(thisMatAllele, split=",")), byrow=T)
  thisPatAllele <- matrix(nrow=length(thisPatAllele), ncol=3, data=unlist(strsplit(thisPatAllele, split=",")), byrow=T)

  thisLoc <- all_significance_effects_d2[i,"POS"]
  if(!thisLoc%in%thisMatAllele[,1]) { problemChildren <- c(problemChildren, i)}
  if(thisLoc%in%thisMatAllele[,1]) {
    matAllele[i] <- thisMatAllele[thisMatAllele[,1]==all_significance_effects_d2[i,"POS"],2]
    patAllele[i] <- thisPatAllele[thisPatAllele[,1]==all_significance_effects_d2[i,"POS"],2]
  }
}

all_significance_effects_d2$MatAllele <- matAllele
all_significance_effects_d2$PatAllele <- patAllele


################################################################################
#Now we add the allele frequency, where data was available.
################################################################################


vep_d2 <- readLines(vep_d2_saveFile)
vep_d2 <- vep_d2[-grep("^##", vep_d2)]
vep_d2_header <- vep_d2[grep("^#", vep_d2)]
vep_d2_header <- gsub("^#", "", vep_d2_header)
vep_d2_header <- strsplit(vep_d2_header, split="\t")[[1]]

vep_d2 <- read.delim(vep_d2_saveFile, header=F, sep="\t", stringsAsFactors=F, comment.char="#")
colnames(vep_d2) <- vep_d2_header

vep_d2 <- vep_d2[,c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
      "Existing_variation", "DISTANCE", "STRAND", "VARIANT_CLASS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE",
      "CANONICAL", "gnomADg_AF", "CLIN_SIG", "PUBMED", "CADD_PHRED", "CADD_RAW", "gnomad3_AF", "GerpRS")]


all_significance_effects_d2$name1 <- paste(all_significance_effects_d2$CHROM, "_", all_significance_effects_d2$POS, "_", all_significance_effects_d2$MatAllele, "/", all_significance_effects_d2$PatAllele, sep="")
all_significance_effects_d2$name2 <- paste(all_significance_effects_d2$CHROM, "_", all_significance_effects_d2$POS, "_", all_significance_effects_d2$PatAllele, "/", all_significance_effects_d2$MatAllele, sep="")


match1 <- match(all_significance_effects_d2$name1, vep_d2$Uploaded_variation)
match2 <- match(all_significance_effects_d2$name2, vep_d2$Uploaded_variation)



Alt_AF1 <- rep(NA, nrow(all_significance_effects_d2))
Alt_AF1 <- vep_d2[match1, "gnomad3_AF"]

Alt_AF2 <- rep(NA, nrow(all_significance_effects_d2))
Alt_AF2 <- vep_d2[match2, "gnomad3_AF"]

Alt_AF <- cbind(Alt_AF1, Alt_AF2)
Alt_AF <- as.data.frame(Alt_AF)
Alt_AF[,1] <- as.numeric(as.character(Alt_AF[,1]))
Alt_AF[,2] <- as.numeric(as.character(Alt_AF[,2]))

final_Alt_AF <- rep(NA, nrow(Alt_AF))
for (i in 1:nrow(Alt_AF)) {
  if(!is.na(Alt_AF[i,1]) && is.na(Alt_AF[i,2])) { final_Alt_AF[i] <- Alt_AF[i,1]}
  if(is.na(Alt_AF[i,1]) && !is.na(Alt_AF[i,2])) { final_Alt_AF[i] <- Alt_AF[i,2]}
  if(!is.na(Alt_AF[i,1]) && !is.na(Alt_AF[i,2])) {
    print(i)
    break()
  }

}


all_significance_effects_d2$Alt_AF <- final_Alt_AF



name1 <- rep(NA, nrow(all_significance_effects_d2))
name2 <- rep(NA, nrow(all_significance_effects_d2))
for (i in 1:nrow(all_significance_effects_d2)) {
  if(i%%1000==0) {print(paste(i, nrow(all_significance_effects_d2)))}
  twoOptions <- all_significance_effects_d2[i,c("MatAllele", "PatAllele")]
  theAncestralAllele <- all_significance_effects_d2[i,"AncAll"]
  theDerivedAllele <- twoOptions[twoOptions!=theAncestralAllele]
  if(length(theDerivedAllele)==1) {
    this_name1 <- paste(all_significance_effects_d2[i,"CHROM"], "_", all_significance_effects_d2[i,"POS"], "_", theAncestralAllele, "/", theDerivedAllele, sep="")
    name1[i] <- this_name1
    this_name2 <- paste(all_significance_effects_d2[i,"CHROM"], "_", all_significance_effects_d2[i,"POS"], "_", theDerivedAllele, "/", theAncestralAllele, sep="")
    name2[i] <- this_name2
  }

}



all_significance_effects_d2$name1 <- name1
all_significance_effects_d2$name2 <- name2


match1 <- match(all_significance_effects_d2$name1, vep_d2$Uploaded_variation)
match2 <- match(all_significance_effects_d2$name2, vep_d2$Uploaded_variation)


DAF1 <- rep(NA, nrow(all_significance_effects_d2))
DAF1 <- vep_d2[match1, "gnomad3_AF"]

DAF2 <- rep(NA, nrow(all_significance_effects_d2))
DAF2 <- vep_d2[match2, "gnomad3_AF"]

DAF <- cbind(DAF1, DAF2)
DAF <- as.data.frame(DAF)
DAF[,1] <- as.numeric(as.character(DAF[,1]))
DAF[,2] <- as.numeric(as.character(DAF[,2]))

#When the derived is the reference, we need to subtract from 1.
DAF[,2] <- 1-DAF[,2]


finalDAF <- rep(NA, nrow(DAF))
for (i in 1:nrow(DAF)) {
  if(!is.na(DAF[i,1]) && is.na(DAF[i,2])) { finalDAF[i] <- DAF[i,1]}
  if(is.na(DAF[i,1]) && !is.na(DAF[i,2])) { finalDAF[i] <- DAF[i,2]}
  if(!is.na(DAF[i,1]) && !is.na(DAF[i,2])) {
    print(i)
    break()
  }

}


all_significance_effects_d2$DAF <- finalDAF

rm(DAF, DAF1, DAF2, finalDAF, matAllele, patAllele, match1, match2, theMatch, theMats, thePats, vep_d2, vep_d2_header)


all_significance_effects_d2 <- all_significance_effects_d2[!is.na(all_significance_effects_d2$MatAllele),]
all_significance_effects_d2 <- all_significance_effects_d2[!is.na(all_significance_effects_d2$PatAllele),]





################################################################################
################################################################################
#Barplots
################################################################################
################################################################################

################################################################################
#Combine the tables.
################################################################################

all_significance_effects_combined <- rbind(all_significance_effects_d1, all_significance_effects_d2)


saveFile <- paste(outDir, "All_significance_effects_bothDonors_withDAF_and_AncestralAlleles_onlyLowestPval.txt", sep="")
write.table(all_significance_effects_combined, saveFile, row.names=F, col.names=T, sep="\t", quote=F)
#all_significance_effects_combined <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

all_significance_effects_combined <- all_significance_effects_combined[!is.na(all_significance_effects_combined$MatAllele),]
all_significance_effects_combined <- all_significance_effects_combined[!is.na(all_significance_effects_combined$PatAllele),]



################################################################################
#Make a figure relating significance to DAF, extending the tail ends a bit.
#Make it a fraction barplot, so we can see "of the with such-and-such DAF,
#here is the distribution of significance..."
################################################################################



finalTable <- data.frame(frequency_set=c(), significance_set=c(), percentage=c())
finalTable_noInput <- data.frame(frequency_set=c(), significance_set=c(), percentage=c())
all_frequencies <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.2, 0.33, 0.5, 0.66, 0.8, 0.9, 0.99, 0.999, 1)
all_significance <- c(0.05, 0.01, 0.001, 1E-5, 1E-10, 1E-15, 1E-20, 0)

theCounts <- c()

mini_table_1 <- all_significance_effects_combined[is.na(all_significance_effects_combined[,"DAF"]),]
mini_table_1 <- mini_table_1[!is.na(mini_table_1$these_pvals),]

mini_table_1_d1 <- mini_table_1[mini_table_1$Donor=="Donor1",]
mini_table_1_d2 <- mini_table_1[mini_table_1$Donor=="Donor2",]
mini_table_1_d1 <- mini_table_1_d1[!mini_table_1_d1$Regions%in%rownames(significant_input_d1),]
mini_table_1_d2 <- mini_table_1_d2[!mini_table_1_d2$Regions%in%rownames(significant_input_d2),]
#mini_table_1_d1 <- mini_table_1_d1[!mini_table_1_d1$Regions%in%rownames(significant_input_d1_stricter),]
#mini_table_1_d2 <- mini_table_1_d2[!mini_table_1_d2$Regions%in%rownames(significant_input_d2_stricter),]
mini_table_1_noInput <- rbind(mini_table_1_d1, mini_table_1_d2)

theCounts <- c(theCounts, nrow(mini_table_1))
for (j in 1:(length(all_significance)-1)) {
  mini_table_2 <- mini_table_1[mini_table_1[,"these_pvals"]<=all_significance[j],]
  mini_table_2 <- mini_table_2[mini_table_2[,"these_pvals"]>all_significance[j+1],]

  thisPerc <- nrow(mini_table_2)/nrow(mini_table_1)
  thisLine <- data.frame(frequency_set=NA, significance_set=all_significance[j], percentage=thisPerc)
  finalTable <- rbind(finalTable, thisLine)

  mini_table_2_noInput <- mini_table_1_noInput[mini_table_1_noInput[,"these_pvals"]<=all_significance[j],]
  mini_table_2_noInput <- mini_table_2_noInput[mini_table_2_noInput[,"these_pvals"]>all_significance[j+1],]

  thisPerc <- nrow(mini_table_2_noInput)/nrow(mini_table_1_noInput)
  thisLine <- data.frame(frequency_set=NA, significance_set=all_significance[j], percentage=thisPerc)
  finalTable_noInput <- rbind(finalTable_noInput, thisLine)

}


mini_table_1 <- all_significance_effects_combined[all_significance_effects_combined[,"DAF"]==0,]
mini_table_1 <- mini_table_1[!is.na(mini_table_1$DAF),]

mini_table_1_d1 <- mini_table_1[mini_table_1$Donor=="Donor1",]
mini_table_1_d2 <- mini_table_1[mini_table_1$Donor=="Donor2",]
mini_table_1_d1 <- mini_table_1_d1[!mini_table_1_d1$Regions%in%rownames(significant_input_d1),]
mini_table_1_d2 <- mini_table_1_d2[!mini_table_1_d2$Regions%in%rownames(significant_input_d2),]
#mini_table_1_d1 <- mini_table_1_d1[!mini_table_1_d1$Regions%in%rownames(significant_input_d1_stricter),]
#mini_table_1_d2 <- mini_table_1_d2[!mini_table_1_d2$Regions%in%rownames(significant_input_d2_stricter),]
mini_table_1_noInput <- rbind(mini_table_1_d1, mini_table_1_d2)

theCounts <- c(theCounts, nrow(mini_table_1))
for (j in 1:(length(all_significance)-1)) {
  mini_table_2 <- mini_table_1[mini_table_1[,"these_pvals"]<=all_significance[j],]
  mini_table_2 <- mini_table_2[mini_table_2[,"these_pvals"]>all_significance[j+1],]
  thisPerc <- nrow(mini_table_2)/nrow(mini_table_1)
  thisLine <- data.frame(frequency_set=0, significance_set=all_significance[j], percentage=thisPerc)
  finalTable <- rbind(finalTable, thisLine)

  mini_table_2_noInput <- mini_table_1_noInput[mini_table_1_noInput[,"these_pvals"]<=all_significance[j],]
  mini_table_2_noInput <- mini_table_2_noInput[mini_table_2_noInput[,"these_pvals"]>all_significance[j+1],]

  thisPerc <- nrow(mini_table_2_noInput)/nrow(mini_table_1_noInput)
  thisLine <- data.frame(frequency_set=0, significance_set=all_significance[j], percentage=thisPerc)
  finalTable_noInput <- rbind(finalTable_noInput, thisLine)
}



for (i in 2:length(all_frequencies)) {
  mini_table_1 <- all_significance_effects_combined[all_significance_effects_combined[,"DAF"]<all_frequencies[i],]
  mini_table_1 <- mini_table_1[mini_table_1[,"DAF"]>=all_frequencies[i-1],]
  mini_table_1 <- mini_table_1[!is.na(mini_table_1[,"DAF"]),]

  mini_table_1_d1 <- mini_table_1[mini_table_1$Donor=="Donor1",]
  mini_table_1_d2 <- mini_table_1[mini_table_1$Donor=="Donor2",]
  mini_table_1_d1 <- mini_table_1_d1[!mini_table_1_d1$Regions%in%rownames(significant_input_d1),]
  mini_table_1_d2 <- mini_table_1_d2[!mini_table_1_d2$Regions%in%rownames(significant_input_d2),]
  #mini_table_1_d1 <- mini_table_1_d1[!mini_table_1_d1$Regions%in%rownames(significant_input_d1_stricter),]
  #mini_table_1_d2 <- mini_table_1_d2[!mini_table_1_d2$Regions%in%rownames(significant_input_d2_stricter),]
  mini_table_1_noInput <- rbind(mini_table_1_d1, mini_table_1_d2)
  mini_table_1_noInput <- mini_table_1_noInput[!is.na(mini_table_1_noInput[,"DAF"]),]

  theCounts <- c(theCounts, nrow(mini_table_1))
  for (j in 1:(length(all_significance)-1)) {
    mini_table_2 <- mini_table_1[mini_table_1[,"these_pvals"]<=all_significance[j],]
    mini_table_2 <- mini_table_2[mini_table_2[,"these_pvals"]>all_significance[j+1],]
    thisPerc <- nrow(mini_table_2)/nrow(mini_table_1)
    thisLine <- data.frame(frequency_set=all_frequencies[i], significance_set=all_significance[j], percentage=thisPerc)
    finalTable <- rbind(finalTable, thisLine)

    mini_table_2_noInput <- mini_table_1_noInput[mini_table_1_noInput[,"these_pvals"]<=all_significance[j],]
    mini_table_2_noInput <- mini_table_2_noInput[mini_table_2_noInput[,"these_pvals"]>all_significance[j+1],]
    thisPerc <- nrow(mini_table_2_noInput)/nrow(mini_table_1_noInput)
    thisLine <- data.frame(frequency_set=all_frequencies[i], significance_set=all_significance[j], percentage=thisPerc)
    finalTable_noInput <- rbind(finalTable_noInput, thisLine)
  }

}

colnames(finalTable) <- c("frequency", "significance", "percentage")
finalTable[,1] <- factor(finalTable[,1], levels=c(NA, all_frequencies))
finalTable[,2] <- factor(finalTable[,2], levels=all_significance[1:(length(all_significance)-1)])
finalTable[,3] <- as.numeric(as.character(finalTable[,3]))

for (i in unique(as.character(finalTable[,1]))) {
  thisSet <- finalTable[as.character(finalTable[,1])==as.character(i),"percentage"]
  thisSet <- thisSet[!is.na(thisSet)]
  print(sum(thisSet))
}

colnames(finalTable_noInput) <- c("frequency", "significance", "percentage")
finalTable_noInput[,1] <- factor(finalTable_noInput[,1], levels=c(NA, all_frequencies))
finalTable_noInput[,2] <- factor(finalTable_noInput[,2], levels=all_significance[1:(length(all_significance)-1)])
finalTable_noInput[,3] <- as.numeric(as.character(finalTable_noInput[,3]))

for (i in unique(as.character(finalTable_noInput[,1]))) {
  thisSet <- finalTable_noInput[as.character(finalTable_noInput[,1])==as.character(i),"percentage"]
  thisSet <- thisSet[!is.na(thisSet)]
  print(sum(thisSet))
}




################################################################################
#Do the same as above, but for noInput cases.
################################################################################


saveFile <- paste(outDir, "Figure_3B.pdf", sep="")
ggplot(finalTable_noInput[!is.na(finalTable_noInput$frequency),], aes(fill=significance, y=percentage, x=frequency)) +
    geom_bar(position="stack", stat="identity") + theme_classic(base_size=25) +
    xlab("Derived Allele Frequency") + ylab("Fraction of Regions") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="bottom", legend.text = element_text(size=15), legend.title=element_text(size=20))
ggsave(saveFile)









#
