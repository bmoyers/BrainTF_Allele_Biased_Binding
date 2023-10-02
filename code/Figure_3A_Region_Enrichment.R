#!/usr/bin/env R
#Figure_3A_Region_Enrichment.R


################################################################################
#This script is used to produce Figure 3A.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges, ggplot2, viridis, and matrixStats packages are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# all_cCREs: path to cCREs, version 4. ENCODE Accession ENCSR800VNX.
# varSeqs_d1: Path to a fasta of all variant sequences for Donor 1, produced
#     by the script Building_fasta_variant_sequences.sh.
# varSeqs_d2: Path to a fasta of all variant sequences for Donor 2, produced
#     by the script Building_fasta_variant_sequences.sh.
# variants_d1: Path to the VCF for donor 1, provided via the following DOI:
#     https://doi.org/10.7303/syn4921369.
# variants_d2: Path to the VCF for donor 2, provided via the following DOI:
#     https://doi.org/10.7303/syn4921369..
# regions_d1: A table compiling p-value calculations across all experiments for
#     donor 1, produced by the script CompilingHaplotypePvalsTable.R
# depths_m_d1: A table compiling the number of reads corresponding to the
#     number of reads for each variant region in each dataset for haplotype 1
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depths_p_d1: As depths_m_d1, but for the alternate haplotype.
# depths_m_d2: As depths_m_d1, but for donor 2.
# depths_p_d2: As depths_m_d2, but for the alternate haplotype.
# allPeaks: A file containing the paths to all available peaks in the dataset.
#     Can be rewritten to be a directroy containing all of the peak sets.
#
################################################################################


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(matrixStats)
library(ggplot2)
library(viridis)
library(GenomicRanges)

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


#The following function takes a genomicRanges object of variants or
#regions of interest, as well as a genomic ranges object for each of
#Downstream, Intronic, Exonic, 3'UTR, 5'UTR, Upstream, Promoter, and Enhancers.
#It then determines which class each variant or region belongs to, following
#the Heirarchy priority:
#Enhancers, Promoter, Upstream, 5’UTR, 3’UTR, Exon, Intron, Downstream, Intergenic
find_genomic_region_types <- function(theRegions, all_cCREs, theStates) {
  #theRegions <- variant_regions
  #set up the genomicRanges of all_cCREs
  all_cCREs <- as.data.frame(all_cCREs)
  colnames(all_cCREs) <- c("chr", "start", "end", "label1", "label2", "cCRE_Type")
  all_cCREs_gr <- makeGRangesFromDataFrame(all_cCREs)

  #Initially label everything as "intergenic"
  theLabels <- rep("None", length(theRegions))

  #Identify Overlaps
  thisOverlap <- as.data.frame(findOverlaps(theRegions, all_cCREs_gr))

  thisOverlap[,2] <- all_cCREs[thisOverlap[,2],"cCRE_Type"]

  for (j in (length(theStates)-1):1) {
    thisSet <- unique(thisOverlap[thisOverlap[,2]==theStates[j],1])
    theLabels[thisSet] <- theStates[j]
  }
  return(theLabels)
}



#The following function creates a summary table of the number of
#elements overlapping with a given cCRE, as well as the fraction
#of all elements overlapped falling into each category.
createSummaryCountTable <- function(labelList, theNames) {
  uniqueLabels <- unique(unlist(labelList))
  summaryTable_counts <- c()
  for (j in 1:length(labelList)) {
    thisVec <- c()
    thisSet <- labelList[[j]]
    for (i in 1:length(uniqueLabels)) { thisVec <- c(thisVec, length(thisSet[thisSet==uniqueLabels[i]]))}
    summaryTable_counts <- rbind(summaryTable_counts, thisVec)
  }
  colnames(summaryTable_counts) <- uniqueLabels
  rownames(summaryTable_counts) <- theNames

  summaryTable_fractions <- summaryTable_counts
  for (j in 1:nrow(summaryTable_fractions)) { summaryTable_fractions[j,] <- summaryTable_fractions[j,]/sum(summaryTable_fractions[j,])}
  return(list(summaryTable_counts, summaryTable_fractions))
}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
all_cCREs <- args[2]
varSeqs_d1 <- args[3]
varSeqs_d2 <- args[4]
variants_d1 <- args[5]
variants_d2 <- args[6]
regions_d1 <- args[7]
depths_m_d1 <- args[8]
depths_p_d1 <- args[9]
regions_d2 <- args[10]
depths_m_d2 <- args[11]
depths_p_d2 <- args[12]
allPeaks <- args[13]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#all_cCREs <- "/cluster/home/bmoyers/ENCODE_Shannon_CellTypeComparison/Datasets_2022November03/GRCh38-cCREs.V4.bed.gz"
#varSeqs_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#varSeqs_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#variants_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants.vcf.table"
#variants_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.table"
#regions_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#depths_m_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depths_p_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#regions_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#depths_m_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depths_p_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#allPeaks <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Motif_ABC_Analyses/path_to_peaks.txt"



theStates <- c("PLS", "pELS", "dELS", "CA-H3K4me3 ", "CA-CTCF", "CA-TF", "TF", "CA", "None")
custom_col <- c("#F51313", "#F59913", "#F5C813", "#EEA4A0", "#43B6F3", "#56D59C", "#6EA38B", "#808080", "#FFFFFF")



################################################################################
#Load in cCREs
################################################################################


all_cCREs <- read.table(all_cCREs, header=F, sep="\t", stringsAsFactors=F)



################################################################################
#Second, load in the variant regions Data.
#Need to restrict this to specifically heterozygous variant regions.
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

toRemove <- grep("chrY", varSeqs_d1[,1])
if(length(toRemove)>0) {
  varSeqs_d1 <- varSeqs_d1[-toRemove,]
}

toRemove <- grep("chrX", varSeqs_d2[,1])
if(length(toRemove)>0) {
  varSeqs_d2 <- varSeqs_d2[-toRemove,]
}

toRemove <- grep("chrY", varSeqs_d2[,1])
if(length(toRemove)>0) {
  varSeqs_d2 <- varSeqs_d2[-toRemove,]
}

varSeqs_d1 <- varSeqs_d1[grep("HET", varSeqs_d1[,2]),]
varSeqs_d2 <- varSeqs_d2[grep("HET", varSeqs_d2[,2]),]

variant_regions_d1 <- matrix(nrow=nrow(varSeqs_d1), ncol=3, data=unlist(strsplit(varSeqs_d1[,1], split=":|-")), byrow=T)
variant_regions_d1 <- as.data.frame(variant_regions_d1)
colnames(variant_regions_d1) <- c("chr", "start", "end")

variant_regions_d2 <- matrix(nrow=nrow(varSeqs_d2), ncol=3, data=unlist(strsplit(varSeqs_d2[,1], split=":|-")), byrow=T)
variant_regions_d2 <- as.data.frame(variant_regions_d2)
colnames(variant_regions_d2) <- c("chr", "start", "end")



variant_regions <- rbind(variant_regions_d1, variant_regions_d2)

variant_regions <- makeGRangesFromDataFrame(variant_regions, keep.extra.columns=F, ignore.strand=T)
variant_regions <- keepStandardChromosomes(variant_regions)


#Alternatively, read in the individual variants, ignoring their haplotypes...

variants_d1 <- read.delim(variants_d1, header=F, sep="\t", stringsAsFactors=F)
variants_d1 <- as.data.frame(variants_d1[,1:2])
variants_d1 <- cbind(variants_d1, variants_d1[,2])
colnames(variants_d1) <- c("chr", "start", "end")


variants_d2 <- read.delim(variants_d2, header=F, sep="\t", stringsAsFactors=F)
variants_d2 <- as.data.frame(variants_d2[,1:2])
variants_d2 <- cbind(variants_d2, variants_d2[,2])
colnames(variants_d2) <- c("chr", "start", "end")


toRemove <- grep("chrX", variants_d1[,1])
if(length(toRemove)>0) {
  variants_d1 <- variants_d1[-toRemove,]
}

toRemove <- grep("chrX", variants_d2[,1])
if(length(toRemove)>0) {
  variants_d2 <- variants_d2[-toRemove,]
}

variants <- rbind(variants_d1, variants_d2)
variants <- makeGRangesFromDataFrame(variants, keep.extra.columns=F, ignore.strand=T)
variants <- keepStandardChromosomes(variants)
variants <- reduce(variants)

################################################################################
#Third, create a genomicRanges object for significant regions, and pick out
#significant variants from it.
#Significance will be based on summing across tissues.
#Need to make a separate set here that removes anything significant in the input.
################################################################################


regions_d1 <- read.table(regions_d1)
depths_m_d1 <- read.table(depths_m_d1)
depths_p_d1 <- read.table(depths_p_d1)
depths_d1 <- depths_m_d1 + depths_p_d1

toRemove <- grep("chrX", regions_d1[,1])
if(length(toRemove)>0) {
  regions_d1 <- regions_d1[-toRemove,]
  depths_m_d1 <- depths_m_d1[-toRemove,]
  depths_p_d1 <- depths_p_d1[-toRemove,]
  depths_d1 <- depths_d1[-toRemove,]
}


#Remove histones and POLR2.


forbiddenExprs <- colnames(regions_d1)[c(grep("H[3|4]K", colnames(regions_d1)), grep("^POL", colnames(regions_d1)))]

regions_d1 <- regions_d1[,!colnames(regions_d1)%in%forbiddenExprs]
depths_m_d1 <- depths_m_d1[,!colnames(depths_m_d1)%in%forbiddenExprs]
depths_p_d1 <- depths_p_d1[,!colnames(depths_p_d1)%in%forbiddenExprs]
depths_d1 <- depths_d1[,!colnames(depths_d1)%in%forbiddenExprs]

#Identify specific cases where the input was significant.

sig_regions_d1_input <- rownames(regions_d1[regions_d1[,"INPUT"]<=0.001,])

#Remove input.


forbiddenExprs <- colnames(regions_d1)[c(grep("INPUT", colnames(regions_d1)))]

regions_d1 <- regions_d1[,!colnames(regions_d1)%in%forbiddenExprs]
depths_m_d1 <- depths_m_d1[,!colnames(depths_m_d1)%in%forbiddenExprs]
depths_p_d1 <- depths_p_d1[,!colnames(depths_p_d1)%in%forbiddenExprs]
depths_d1 <- depths_d1[,!colnames(depths_d1)%in%forbiddenExprs]


#Identify all significant regions.
#Remove cases where input was also significant.

sig_regions_d1 <- rownames(regions_d1)[rowMins(as.matrix(regions_d1))<0.001]
sig_regions_d1_noInput <- sig_regions_d1[!sig_regions_d1%in%sig_regions_d1_input]
sig_regions_d1 <- matrix(nrow=length(sig_regions_d1), ncol=3, data=unlist(strsplit(sig_regions_d1, split=":|-")), byrow=T)
sig_regions_d1 <- as.data.frame(sig_regions_d1)
colnames(sig_regions_d1) <- c("chr", "start", "end")

sig_regions_d1_noInput <- matrix(nrow=length(sig_regions_d1_noInput), ncol=3, data=unlist(strsplit(sig_regions_d1_noInput, split=":|-")), byrow=T)
sig_regions_d1_noInput <- as.data.frame(sig_regions_d1_noInput)
colnames(sig_regions_d1_noInput) <- c("chr", "start", "end")


#Identify all regions with sufficient depth.
possible_regions_d1 <- rownames(depths_d1)[rowMaxs(as.matrix(depths_d1))>=11]
possible_regions_d1 <- matrix(nrow=length(possible_regions_d1), ncol=3, data=unlist(strsplit(possible_regions_d1, split=":|-")), byrow=T)
possible_regions_d1 <- as.data.frame(possible_regions_d1)
colnames(possible_regions_d1) <- c("chr", "start", "end")


################################################################################
#All the same things, but for donor 2.
################################################################################


regions_d2 <- read.table(regions_d2)
depths_m_d2 <- read.table(depths_m_d2)
depths_p_d2 <- read.table(depths_p_d2)
depths_d2 <- depths_m_d2 + depths_p_d2


toRemove <- grep("chrX", regions_d2[,1])
if(length(toRemove)>0) {
  regions_d2 <- regions_d2[-toRemove,]
  depths_m_d2 <- depths_m_d2[-toRemove,]
  depths_p_d2 <- depths_p_d2[-toRemove,]
  depths_d2 <- depths_d2[-toRemove,]
}

#Remove histones and POLR2.


forbiddenExprs <- colnames(regions_d2)[c(grep("H[3|4]K", colnames(regions_d2)), grep("^POL", colnames(regions_d2)))]

regions_d2 <- regions_d2[,!colnames(regions_d2)%in%forbiddenExprs]
depths_m_d2 <- depths_m_d2[,!colnames(depths_m_d2)%in%forbiddenExprs]
depths_p_d2 <- depths_p_d2[,!colnames(depths_p_d2)%in%forbiddenExprs]
depths_d2 <- depths_d2[,!colnames(depths_d2)%in%forbiddenExprs]


#Identify specific cases where the input was significant.

sig_regions_d2_input <- rownames(regions_d2[regions_d2[,"INPUT"]<=0.001,])

#Remove input.


forbiddenExprs <- colnames(regions_d2)[c(grep("INPUT", colnames(regions_d2)))]

regions_d2 <- regions_d2[,!colnames(regions_d2)%in%forbiddenExprs]
depths_m_d2 <- depths_m_d2[,!colnames(depths_m_d2)%in%forbiddenExprs]
depths_p_d2 <- depths_p_d2[,!colnames(depths_p_d2)%in%forbiddenExprs]
depths_d2 <- depths_d2[,!colnames(depths_d2)%in%forbiddenExprs]


#Identify all significant regions.
#Remove cases where input was also significant.


sig_regions_d2 <- rownames(regions_d2)[rowMins(as.matrix(regions_d2))<0.001]
sig_regions_d2_noInput <- sig_regions_d2[!sig_regions_d2%in%sig_regions_d2_input]
sig_regions_d2 <- matrix(nrow=length(sig_regions_d2), ncol=3, data=unlist(strsplit(sig_regions_d2, split=":|-")), byrow=T)
sig_regions_d2 <- as.data.frame(sig_regions_d2)
colnames(sig_regions_d2) <- c("chr", "start", "end")


sig_regions_d2_noInput <- matrix(nrow=length(sig_regions_d2_noInput), ncol=3, data=unlist(strsplit(sig_regions_d2_noInput, split=":|-")), byrow=T)
sig_regions_d2_noInput <- as.data.frame(sig_regions_d2_noInput)
colnames(sig_regions_d2_noInput) <- c("chr", "start", "end")


possible_regions_d2 <- rownames(depths_d2)[rowMaxs(as.matrix(depths_d2))>=11]
possible_regions_d2 <- matrix(nrow=length(possible_regions_d2), ncol=3, data=unlist(strsplit(possible_regions_d2, split=":|-")), byrow=T)
possible_regions_d2 <- as.data.frame(possible_regions_d2)
colnames(possible_regions_d2) <- c("chr", "start", "end")


################################################################################
#Combine the two donors.
################################################################################

sig_regions <- rbind(sig_regions_d1, sig_regions_d2)
sig_regions <- makeGRangesFromDataFrame(sig_regions, ignore.strand=T, keep.extra.columns=F)
sig_regions <- keepStandardChromosomes(sig_regions)

sig_regions_noInput <- rbind(sig_regions_d1_noInput, sig_regions_d2_noInput)
sig_regions_noInput <- makeGRangesFromDataFrame(sig_regions_noInput, ignore.strand=T, keep.extra.columns=F)
sig_regions_noInput <- keepStandardChromosomes(sig_regions_noInput)

possible_regions <- rbind(possible_regions_d1, possible_regions_d1)
possible_regions <- makeGRangesFromDataFrame(possible_regions, ignore.strand=T, keep.extra.columns=F)
possible_regions <- keepStandardChromosomes(possible_regions)

theIntersect <- as.matrix(findOverlaps(variants, sig_regions))
sig_variants <- as.data.frame(variants)[unique(theIntersect[,1]),]
sig_variants <- makeGRangesFromDataFrame(sig_variants, ignore.strand=T, keep.extra.columns=F)
sig_variants <- keepStandardChromosomes(sig_variants)
sig_variants <- reduce(sig_variants)

theIntersect <- as.matrix(findOverlaps(variants, sig_regions_noInput))
sig_variants_noInput <- as.data.frame(variants)[unique(theIntersect[,1]),]
sig_variants_noInput <- makeGRangesFromDataFrame(sig_variants_noInput, ignore.strand=T, keep.extra.columns=F)
sig_variants_noInput <- keepStandardChromosomes(sig_variants_noInput)
sig_variants_noInput <- reduce(sig_variants_noInput)

theIntersect <- as.matrix(findOverlaps(variants, possible_regions))
possible_variants <- as.data.frame(variants)[unique(theIntersect[,1]),]
possible_variants <- makeGRangesFromDataFrame(possible_variants, ignore.strand=T, keep.extra.columns=F)
possible_variants <- keepStandardChromosomes(possible_variants)


################################################################################
#Fourth, we're going to find all the peaks for all experiments and check out
#their locations and distributions...
################################################################################

allPeaks <- readLines(allPeaks)

################################################################################
#For each file identify the TF and the tissue to make sure
#everything is properly formatted and we don't have weird outliers...
################################################################################

theTFs <- rep(NA, length(allPeaks))
theTissues <-  rep(NA, length(allPeaks))
for (i in 1:length(allPeaks)) {
  thisSet <- strsplit(allPeaks[i], split="peak/spp/idr/optimal_set/", fixed=T)[[1]][2]
  thisSet <- strsplit(thisSet, split="_")[[1]]
  theTFs[i] <- thisSet[1]
  theTissues[i] <- thisSet[2]
}
theTissues <- toupper(theTissues)

allPeaks <- cbind(theTFs, theTissues, allPeaks)

allPeaks <- allPeaks[!allPeaks[,"theTFs"]%in%c("H3K27ac", "H3K4me3", "H3K9ac", "Pol2"),]


full_peakSet <- c()
for (i in 1:nrow(allPeaks)) {
  thesePeaks <- read.table(allPeaks[i,3], header=F, sep="\t", stringsAsFactors=F)
  thesePeaks$TF <- rep(allPeaks[i,1], nrow(thesePeaks))
  thesePeaks$Tissue <- rep(allPeaks[i,2], nrow(thesePeaks))
  theCenters <- thesePeaks[,2]+thesePeaks[,10]-1
  theStarts <- theCenters-50
  theStops <- theCenters+50
  thesePeaks[,2] <- theStarts
  thesePeaks[,3] <- theStops
  full_peakSet <- rbind(full_peakSet, thesePeaks)
}



colnames(full_peakSet)[1:3] <- c("chr", "start", "end")
allPeaks <- makeGRangesFromDataFrame(full_peakSet)
allPeaks <- keepStandardChromosomes(allPeaks)

toRemove <- as.data.frame(findOverlaps(allPeaks, sig_regions))

allPeaks_noSig <- allPeaks[-unique(toRemove[,1])]
allPeaks_noSig <- reduce(allPeaks_noSig)
allPeaks_noSig <- as.data.frame(allPeaks_noSig)
saveFile <- paste(outDir, "Background_regions_for_GREAT_analysis.bed", sep="")
write.table(allPeaks_noSig[,1:3], saveFile, row.names=F, col.names=F, sep="\t", quote=F)


allPeaks <- reduce(allPeaks)


################################################################################
#Fifth, we're going to subset the variant (and sig variant) regions into those
#found in peaks.
################################################################################

theIntersect <- as.matrix(findOverlaps(variant_regions, allPeaks))
peak_variant_regions <- as.data.frame(variant_regions)[unique(theIntersect[,1]),]
peak_variant_regions <- peak_variant_regions[,1:3]
peak_variant_regions <- makeGRangesFromDataFrame(peak_variant_regions, keep.extra.columns=F, ignore.strand=T)
peak_variant_regions <- keepStandardChromosomes(peak_variant_regions)

theIntersect <- as.matrix(findOverlaps(variants, allPeaks))
peak_variants <- as.data.frame(variants)[unique(theIntersect[,1]),]
peak_variants <- peak_variants[,1:3]
peak_variants <- makeGRangesFromDataFrame(peak_variants, keep.extra.columns=F, ignore.strand=T)
peak_variants <- keepStandardChromosomes(peak_variants)

theIntersect <- as.matrix(findOverlaps(sig_regions, allPeaks))
peak_sig_regions <- as.data.frame(sig_regions)[unique(theIntersect[,1]),]
peak_sig_regions <- peak_sig_regions[,1:3]
peak_sig_regions <- makeGRangesFromDataFrame(peak_sig_regions, keep.extra.columns=F, ignore.strand=T)
peak_sig_regions <- keepStandardChromosomes(peak_sig_regions)

theIntersect <- as.matrix(findOverlaps(sig_regions_noInput, allPeaks))
peak_sig_regions_noInput <- as.data.frame(sig_regions_noInput)[unique(theIntersect[,1]),]
peak_sig_regions_noInput <- peak_sig_regions_noInput[,1:3]
peak_sig_regions_noInput <- makeGRangesFromDataFrame(peak_sig_regions_noInput, keep.extra.columns=F, ignore.strand=T)
peak_sig_regions_noInput <- keepStandardChromosomes(peak_sig_regions_noInput)

theIntersect <- as.matrix(findOverlaps(sig_variants, allPeaks))
peak_sig_variants <- as.data.frame(sig_variants)[unique(theIntersect[,1]),]
peak_sig_variants <- peak_sig_variants[,1:3]
peak_sig_variants <- makeGRangesFromDataFrame(peak_sig_variants, keep.extra.columns=F, ignore.strand=T)
peak_sig_variants <- keepStandardChromosomes(peak_sig_variants)
peak_sig_variants <- reduce(peak_sig_variants)

theIntersect <- as.matrix(findOverlaps(sig_variants_noInput, allPeaks))
peak_sig_variants_noInput <- as.data.frame(sig_variants_noInput)[unique(theIntersect[,1]),]
peak_sig_variants_noInput <- peak_sig_variants_noInput[,1:3]
peak_sig_variants_noInput <- makeGRangesFromDataFrame(peak_sig_variants_noInput, keep.extra.columns=F, ignore.strand=T)
peak_sig_variants_noInput <- keepStandardChromosomes(peak_sig_variants_noInput)
peak_sig_variants_noInput <- reduce(peak_sig_variants_noInput)



################################################################################
#Side-task, for all significant variants, determine the percentage overlap with
#cCREs.
################################################################################
colnames(all_cCREs)[1:3] <- c("chr", "start", "end")
all_cCREs_gr <- makeGRangesFromDataFrame(all_cCREs)
sidetaskIntersect <- as.data.frame(findOverlaps(sig_regions_noInput, all_cCREs_gr))
length(unique(sidetaskIntersect[,1]))
1 - length(unique(sidetaskIntersect[,1])) / length(sig_regions_noInput)



################################################################################
#Side-task 2; For each TF, determine the regions in which it is significant,
#and which of those has variants within a peak. For each donor.
################################################################################

uniqueTFs <- unique(full_peakSet$TF)

fraction_sig_in_peaks_table <- c()
for (i in 1:length(uniqueTFs)) {
  thesePeaks <- full_peakSet[full_peakSet$TF==uniqueTFs[i],]
  thesePeaks_gr <- makeGRangesFromDataFrame(thesePeaks)

  this_totalSig <- 0
  this_sigInPeaks <- 0
  if(uniqueTFs[i]%in%colnames(regions_d1)) {
    these_sig_d1 <- rownames(regions_d1)[which(regions_d1[,uniqueTFs[i]]<=0.001)]
    these_sig_d1 <- as.data.frame(matrix(nrow=length(these_sig_d1), ncol=3, data=unlist(strsplit(these_sig_d1, split=":|-")), byrow=T))
    colnames(these_sig_d1) <- c("chr", "start", "end")
    these_sig_d1[,2] <- as.numeric(as.character(these_sig_d1[,2]))
    these_sig_d1[,3] <- as.numeric(as.character(these_sig_d1[,3]))
    these_sig_d1[,2] <- these_sig_d1[,2]+100
    these_sig_d1[,3] <- these_sig_d1[,3]-100
    these_sig_d1_gr <- makeGRangesFromDataFrame(these_sig_d1)
    this_sigPeaks_intersect <- as.data.frame(findOverlaps(these_sig_d1_gr, thesePeaks_gr))
    this_totalSig <- this_totalSig + length(these_sig_d1_gr)
    this_sigInPeaks <- this_sigInPeaks + length(unique(this_sigPeaks_intersect[,1]))
  }
  if(uniqueTFs[i]%in%colnames(regions_d2)) {
    these_sig_d2 <- rownames(regions_d2)[which(regions_d2[,uniqueTFs[i]]<=0.001)]
    these_sig_d2 <- as.data.frame(matrix(nrow=length(these_sig_d2), ncol=3, data=unlist(strsplit(these_sig_d2, split=":|-")), byrow=T))
    colnames(these_sig_d2) <- c("chr", "start", "end")
    these_sig_d2[,2] <- as.numeric(as.character(these_sig_d2[,2]))
    these_sig_d2[,3] <- as.numeric(as.character(these_sig_d2[,3]))
    these_sig_d2[,2] <- these_sig_d2[,2]+100
    these_sig_d2[,3] <- these_sig_d2[,3]-100
    these_sig_d2_gr <- makeGRangesFromDataFrame(these_sig_d2)
    this_sigPeaks_intersect <- as.data.frame(findOverlaps(these_sig_d2_gr, thesePeaks_gr))
    this_totalSig <- this_totalSig + length(these_sig_d2_gr)
    this_sigInPeaks <- this_sigInPeaks + length(unique(this_sigPeaks_intersect[,1]))
  }

  thisFraction <- this_sigInPeaks / this_totalSig
  fraction_sig_in_peaks_table <- rbind(fraction_sig_in_peaks_table, c(uniqueTFs[i], this_sigInPeaks, this_totalSig, thisFraction))
}

fraction_sig_in_peaks_table <- as.data.frame(fraction_sig_in_peaks_table)
colnames(fraction_sig_in_peaks_table) <- c("TF", "Sig_in_Peaks", "totalSig", "Fraction")
for (i in 2:ncol(fraction_sig_in_peaks_table)) { fraction_sig_in_peaks_table[,i] <- as.numeric(as.character(fraction_sig_in_peaks_table[,i]))}

1 - sum(fraction_sig_in_peaks_table$Sig_in_Peaks) / sum(fraction_sig_in_peaks_table$totalSig)


################################################################################
#Sixth, perform heirarchical overlaps.
################################################################################

variantRegion_labels <- find_genomic_region_types(variant_regions, all_cCREs, theStates)
variant_labels <- find_genomic_region_types(variants, all_cCREs, theStates)
sigRegions_labels <- find_genomic_region_types(sig_regions, all_cCREs, theStates)
sigRegions_noInput_labels <- find_genomic_region_types(sig_regions_noInput, all_cCREs, theStates)

sigVariants_labels <- find_genomic_region_types(sig_variants, all_cCREs, theStates)
sigVariants_noInput_labels <- find_genomic_region_types(sig_variants_noInput, all_cCREs, theStates)

possibleRegions_labels <- find_genomic_region_types(possible_regions, all_cCREs, theStates)
possibleVariants_labels <- find_genomic_region_types(possible_variants, all_cCREs, theStates)
peakRegions_labels <- find_genomic_region_types(peak_variant_regions, all_cCREs, theStates)

peakVariants_labels <- find_genomic_region_types(peak_variants, all_cCREs, theStates)
peakSigRegions_labels <- find_genomic_region_types(peak_sig_regions, all_cCREs, theStates)
peakSigRegions_noInput_labels <- find_genomic_region_types(peak_sig_regions_noInput, all_cCREs, theStates)
peakSigVariants_labels <- find_genomic_region_types(peak_sig_variants, all_cCREs, theStates)
peakSigVariants_noInput_labels <- find_genomic_region_types(peak_sig_variants_noInput, all_cCREs, theStates)

allPeaks_labels <- find_genomic_region_types(allPeaks, all_cCREs, theStates)



labelList <- list(variantRegion_labels, variant_labels, sigRegions_labels, sigRegions_noInput_labels,
    sigVariants_labels, sigVariants_noInput_labels, possibleRegions_labels, possibleVariants_labels,
    peakRegions_labels, peakVariants_labels, peakSigRegions_labels, peakSigRegions_noInput_labels,
    peakSigVariants_labels, peakSigVariants_noInput_labels, allPeaks_labels)
#

theNames <- list("variantRegion", "variant", "sigRegions", "sigRegions_noInput", "sigVariants",
    "sigVariants_noInput", "possibleRegions", "possibleVariants", "peakRegions", "peakVariants",
    "peakSigRegions", "peakSigRegions_noInput", "peakSigVariants", "peakSigVariants_noInput",
    "allPeaks")

tableList <-  createSummaryCountTable(labelList, theNames)
summaryTable_counts <- tableList[[1]]
summaryTable_fractions <- tableList[[2]]

saveFile <- paste(outDir, "cCRE_Type_Distribution_table_counts.txt", sep="")
write.table(summaryTable_counts, saveFile, row.names=T, col.names=T, sep="\t", quote=F)

saveFile <- paste(outDir, "cCRE_Distribution_table_fractions.txt", sep="")
write.table(summaryTable_fractions, saveFile, row.names=T, col.names=T, sep="\t", quote=F)




################################################################################
#Chi-sq test, of sigPeaks versus peaks.
################################################################################
a <- chisq.test(t(as.matrix(summaryTable_counts[c("peakSigRegions", "peakRegions"),])))
#Highly significant, of course.

################################################################################
#Can also do a binomial test on the PLS, specifically.
################################################################################

pls_prob <- summaryTable_fractions["peakRegions", "PLS"]
pls_numObs <- summaryTable_counts["peakSigRegions","PLS"]
totalNum <- sum(summaryTable_counts["peakSigRegions",])
b <- binom.test(x=pls_numObs, n=totalNum, p=pls_prob)
#p-value < 2.2e-16
#0.1928293 versus 0.1042856
#ratio: 1.84905

pels_prob <- summaryTable_fractions["peakRegions", "pELS"]
pels_numObs <- summaryTable_counts["peakSigRegions","pELS"]
totalNum <- sum(summaryTable_counts["peakSigRegions",])
b <- binom.test(x=pels_numObs, n=totalNum, p=pels_prob)
#p-value = 9.135e-08
#0.1391001 versus 0.1585876
#ratio: 0.8771184

dels_prob <- summaryTable_fractions["peakRegions", "dELS"]
dels_numObs <- summaryTable_counts["peakSigRegions","dELS"]
totalNum <- sum(summaryTable_counts["peakSigRegions",])
b <- binom.test(x=dels_numObs, n=totalNum, p=dels_prob)






################################################################################
#Create a graph of the relevant results.
################################################################################

allRegionsList <- list(peakSigRegions_labels, sigRegions_labels, allPeaks_labels, variantRegion_labels)
names(allRegionsList) <- c("Signif. In Peaks", "Signif.", "Peaks", "Haplotype")

barplot_df <- c()
for (i in 1:length(allRegionsList)) {
  thisSet <- table(allRegionsList[[i]])/length(allRegionsList[[i]])
  for (j in 1:length(thisSet)) {
    barplot_df <- rbind(barplot_df, c(thisSet[j], names(thisSet)[j], names(allRegionsList)[i]))
  }

}
barplot_df <- as.data.frame(barplot_df)
barplot_df[,1] <- as.numeric(as.character(barplot_df[,1]))
barplot_df[,2] <- as.character(barplot_df[,2])
barplot_df[,3] <- as.character(barplot_df[,3])
colnames(barplot_df) <- c("Fraction", "Region", "Dataset")
barplot_df[,3] <- factor(barplot_df[,3], levels=c("Signif. In Peaks", "Signif.", "Peaks", "Haplotype"))
barplot_df[,2] <- factor(barplot_df[,2], levels=theStates)




#saveFile <- paste(outDir, "cCRE_Region_Enrichment_Barplot.pdf", sep="")
#ggplot(barplot_df, aes(fill=Region, y=Fraction, x=Dataset)) +
#    geom_bar(position="fill", stat="identity") + coord_flip() + theme_classic(base_size = 25) +
#    scale_fill_manual(name="cCREs", values = custom_col, na.value="grey50") + xlab("") +
#    theme(legend.position="bottom", legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12))
#ggsave(saveFile, width=8)






################################################################################
#Same, but after removing cases where the Input was significant.
################################################################################

allRegionsList <- list(peakSigRegions_noInput_labels, sigRegions_noInput_labels, allPeaks_labels, variantRegion_labels)
names(allRegionsList) <- c("Signif. In Peaks", "Signif.", "Peaks", "Haplotype")

barplot_df <- c()
for (i in 1:length(allRegionsList)) {
  thisSet <- table(allRegionsList[[i]])/length(allRegionsList[[i]])
  for (j in 1:length(thisSet)) {
    barplot_df <- rbind(barplot_df, c(thisSet[j], names(thisSet)[j], names(allRegionsList)[i]))
  }

}
barplot_df <- as.data.frame(barplot_df)
barplot_df[,1] <- as.numeric(as.character(barplot_df[,1]))
barplot_df[,2] <- as.character(barplot_df[,2])
barplot_df[,3] <- as.character(barplot_df[,3])
colnames(barplot_df) <- c("Fraction", "Region", "Dataset")
barplot_df[,3] <- factor(barplot_df[,3], levels=c("Signif. In Peaks", "Signif.", "Peaks", "Haplotype"))
barplot_df[,2] <- factor(barplot_df[,2], levels=theStates)



saveFile <- paste(outDir, "Figure_3A.pdf", sep="")
ggplot(barplot_df, aes(fill=Region, y=Fraction, x=Dataset)) +
    geom_bar(position="fill", stat="identity") + coord_flip() + theme_classic(base_size = 25) +
    scale_fill_manual(name="cCREs", values = custom_col, na.value="grey50") + xlab("") +
    theme(legend.position="bottom", legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12))
ggsave(saveFile, width=8)







#
