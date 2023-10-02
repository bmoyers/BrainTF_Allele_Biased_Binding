#!/usr/bin/env R
#Figure_4B_4C_4D_5.R



################################################################################
#This script is used to produce Figures 4B-F and Supplemental Figure 11.
#
#Run under R 4.1.0, but can be used under other versions so long as the
#     GenomicRanges, ggplot2, and matrixStats packages are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# d1_fullGTEX_inOurSet_saveFile: Table with phase-linked RNA and TF biased
#     variants for donor 1, produced by the script Finding_significant_Expression_bias_and_linked_TF_bias.R
# d2_fullGTEX_inOurSet_saveFile: Table with phase-linked RNA and TF biased
#     variants for donor 2, produced by the script Finding_significant_Expression_bias_and_linked_TF_bias.R
# fimo_d1_mat: fimo results for haplotype 1, donor 1. Produced by the script
#     Figure_3D_Supplemental_10.R
# fimo_d1_pat: fimo results for haplotype 2, donor 1. Produced by the script
#     Figure_3D_Supplemental_10.R
# fimo_d2_mat: fimo results for haplotype 1, donor 2. Produced by the script
#     Figure_3D_Supplemental_10.R
# fimo_d2_pat: fimo results for haplotype 2, donor 2. Produced by the script
#     Figure_3D_Supplemental_10.R
#
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
library(gridExtra)
library(memes)
library(universalmotif)
library(plotgardener)
library(EnsDb.Hsapiens.v86)
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")


################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


################################################################################
#The following function determines changes in motif strength between alleles.
#1) First, for a given variant sequence, identify all of the motifs which were found
#   in said sequence in both the maternal and paternal.
#2) Second, identify those which were entirely removed or gained. Save in a table.
#3) Third, identify cases of the same motif which has had a change in score, and save
#   in a table.
################################################################################
determine_motif_destruction <- function(maternal_Hits, paternal_Hits, relevantSeqs_gr, relevantSeqs) {
  #maternal_Hits <- fimo_d1_mat
  #paternal_Hits <- fimo_d1_pat
  maternal_Hits_gr <- makeGRangesFromDataFrame(maternal_Hits)
  paternal_Hits_gr <- makeGRangesFromDataFrame(paternal_Hits)

  theIntersect_maternal <- as.data.frame(findOverlaps(maternal_Hits_gr, relevantSeqs_gr))
  theIntersect_paternal <- as.data.frame(findOverlaps(paternal_Hits_gr, relevantSeqs_gr))

  uniqueRegionIndices <- unique(c(unique(theIntersect_maternal[,2]), unique(theIntersect_paternal[,2])))
  uniqueRegionIndices <- uniqueRegionIndices[order(uniqueRegionIndices)]

  motif_destruction_info_table <- c()
  for (i in 1:length(uniqueRegionIndices)) {
    thisSet_maternal <- as.data.frame(maternal_Hits[theIntersect_maternal[theIntersect_maternal[,2]==uniqueRegionIndices[i],1],])
    thisSet_paternal <- as.data.frame(paternal_Hits[theIntersect_paternal[theIntersect_paternal[,2]==uniqueRegionIndices[i],1],])

    unique_maternal_names <- paste(thisSet_maternal[,1], thisSet_maternal[,4], thisSet_maternal[,5], thisSet_maternal[,6], thisSet_maternal[,7], thisSet_maternal[,8], thisSet_maternal[,9], thisSet_maternal[,10], thisSet_maternal[,11], sep="_")
    unique_paternal_names <- paste(thisSet_paternal[,1], thisSet_paternal[,4], thisSet_paternal[,5], thisSet_paternal[,6], thisSet_paternal[,7], thisSet_paternal[,8], thisSet_paternal[,9], thisSet_paternal[,10], thisSet_paternal[,11], sep="_")

    ##############################################################################
    #Remove cases where the motif is exactly the same score, regardless of location.
    #note that if we have repetitive sequence, this might cause some issues.
    #But we want to deemphasize those sequences for this analysis anyway.
    #This removes cases where indels make the appearance of novel motifs when
    #in reality they are the same motif.
    ##############################################################################
    thisSet_maternal <- thisSet_maternal[which(!unique_maternal_names%in%unique_paternal_names),]
    thisSet_paternal <- thisSet_paternal[which(!unique_paternal_names%in%unique_maternal_names),]

    if(nrow(thisSet_maternal)>0 || nrow(thisSet_paternal)>0) {
      #print(i)
      #break()
      if(nrow(thisSet_maternal)>0 && nrow(thisSet_paternal)>0) {
        ##############################################################################
        #I now have all cases in which there has been a change.
        #from here I want to identify whether the change was
        #1) A complete destruction/creation of the motif or
        #2) A weakinging/strengthening of the motif.
        #I can do this by determining if the location and motifID are the same.
        #   between any two entries. If they are, it is a weakening/strengthening.
        #   If there is a case where an entry is not matched between the two sets,
        #   it is a destruction/creation.
        ##############################################################################

        ##############################################################################
        #At this point, we re-introduce location information so we can compare the scores
        #of motifs at the same location when the motif itself has had a change in its sequence.
        #Because we want to compare scores, though, we leave that information off of these names.
        ##############################################################################

        unique_maternal_names <- paste(thisSet_maternal[,1], thisSet_maternal[,2], thisSet_maternal[,3], thisSet_maternal[,4], thisSet_maternal[,5], thisSet_maternal[,6], sep="_")
        unique_paternal_names <- paste(thisSet_paternal[,1], thisSet_paternal[,2], thisSet_paternal[,3], thisSet_paternal[,4], thisSet_paternal[,5], thisSet_paternal[,6], sep="_")

        ##############################################################################
        #If there was simply a weakening/strengthening
        ##############################################################################
        common_names <- unique_maternal_names[unique_maternal_names%in%unique_paternal_names]
        if(length(common_names)>0) {
          for (j in 1:length(common_names)) {
            this_maternal_info <- thisSet_maternal[unique_maternal_names==common_names[j],1:9]
            this_paternal_info <- thisSet_paternal[unique_paternal_names==common_names[j],8:9]
            colnames(this_maternal_info)[8:9] <- paste(colnames(this_maternal_info)[8:9], "maternal", sep="_")
            colnames(this_paternal_info) <- paste(colnames(this_paternal_info), "paternal", sep="_")
            thisLine <- cbind(this_maternal_info, this_paternal_info)

            thisLine_gr <- makeGRangesFromDataFrame(thisLine)
            lineIntersect <- as.data.frame(findOverlaps(thisLine_gr,relevantSeqs_gr))
            thisLine$region <- paste(relevantSeqs[lineIntersect[1,2],1], ":", relevantSeqs[lineIntersect[1,2],2], "-", relevantSeqs[lineIntersect[1,2],3], sep="")

            motif_destruction_info_table <- rbind(motif_destruction_info_table, thisLine)

          }
          thisSet_maternal <- thisSet_maternal[!unique_maternal_names%in%common_names,]
          thisSet_paternal <- thisSet_paternal[!unique_paternal_names%in%common_names,]
        }
        ##############################################################################
        #If there was a creation/destruction
        ##############################################################################

        if(nrow(thisSet_maternal)>0) {
          for (j in 1:nrow(thisSet_maternal)) {
            thisLine <- thisSet_maternal[j,1:9]
            colnames(thisLine)[8:9] <- paste(colnames(thisLine)[8:9], "_maternal", sep="")
            thisLine$score_paternal <- NA
            thisLine$pvalue_paternal <- NA

            thisLine_gr <- makeGRangesFromDataFrame(thisLine)
            lineIntersect <- as.data.frame(findOverlaps(thisLine_gr,relevantSeqs_gr))
            thisLine$region <- paste(relevantSeqs[lineIntersect[1,2],1], ":", relevantSeqs[lineIntersect[1,2],2], "-", relevantSeqs[lineIntersect[1,2],3], sep="")

            motif_destruction_info_table <- rbind(motif_destruction_info_table, thisLine)
          }
        }
        if(nrow(thisSet_paternal)>0) {
          for (j in 1:nrow(thisSet_paternal)) {
            thisLine <- thisSet_paternal[j,1:7]
            thisLine$score_maternal <- NA
            thisLine$pvalue_maternal <- NA
            thisLine$score_paternal <- thisSet_paternal[j,8]
            thisLine$pvalue_paternal <- thisSet_paternal[j,9]

            thisLine_gr <- makeGRangesFromDataFrame(thisLine)
            lineIntersect <- as.data.frame(findOverlaps(thisLine_gr,relevantSeqs_gr))
            thisLine$region <- paste(relevantSeqs[lineIntersect[1,2],1], ":", relevantSeqs[lineIntersect[1,2],2], "-", relevantSeqs[lineIntersect[1,2],3], sep="")

            motif_destruction_info_table <- rbind(motif_destruction_info_table, thisLine)
          }
        }
      }

    }
  }
  return(motif_destruction_info_table)

}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
d1_fullGTEX_inOurSet_saveFile <- args[2]
d2_fullGTEX_inOurSet_saveFile <- args[3]
fimo_d1_mat <- args[4]
fimo_d1_pat <- args[5]
fimo_d2_mat <- args[6]
fimo_d2_pat <- args[7]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#d1_fullGTEX_inOurSet_saveFile <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/GTEX_with_rnaTFInfo_donor1_brainOnly.txt"
#d2_fullGTEX_inOurSet_saveFile <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/GTEX_with_rnaTFInfo_donor2_brainOnly.txt"
#fimo_d1_mat <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/JASPAR_fimo_donor1_maternal.txt"
#fimo_d1_pat <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/JASPAR_fimo_donor1_paternal.txt"
#fimo_d2_mat <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/JASPAR_fimo_donor2_maternal.txt"
#fimo_d2_pat <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/JASPAR_fimo_donor2_paternal.txt"



################################################################################
#Load in linked TF bias and RNA bias files. Create Genomic Ranges
################################################################################

d1_fullGTEX_inOurSet <- read.table(d1_fullGTEX_inOurSet_saveFile, header=T, sep="\t", stringsAsFactors=F)
d2_fullGTEX_inOurSet <- read.table(d2_fullGTEX_inOurSet_saveFile, header=T, sep="\t", stringsAsFactors=F)


##########################################################################
#Restrict to cases within 2kb of the relevant tss
################################################################################

d1_fullGTEX_inOurSet_candidates <- d1_fullGTEX_inOurSet[abs(d1_fullGTEX_inOurSet$tss_distance)<=2000,]
d2_fullGTEX_inOurSet_candidates <- d2_fullGTEX_inOurSet[abs(d2_fullGTEX_inOurSet$tss_distance)<=2000,]


################################################################################
#Add motif disruption information. This becomes slightly more complicated,
#only because there may in fact be multiple TFs with disrupted motifs for
#a given region.
################################################################################

fimo_d1_mat <- read.table(fimo_d1_mat, header=T, sep="\t", stringsAsFactors=F)
fimo_d1_pat <- read.table(fimo_d1_pat, header=T, sep="\t", stringsAsFactors=F)

relevantSeqs <- unique(matrix(nrow=nrow(d1_fullGTEX_inOurSet_candidates), ncol=3, data=unlist(strsplit(d1_fullGTEX_inOurSet_candidates$tf_region, split=":|-")), byrow=T))
relevantSeqs <- as.data.frame(relevantSeqs)
colnames(relevantSeqs) <- c("chr", "start", "end")
relevantSeqs_gr <- makeGRangesFromDataFrame(relevantSeqs)
relevantSeqs_gr <- unique(relevantSeqs_gr)

motif_destruction_info_table_d1 <- determine_motif_destruction(fimo_d1_mat, fimo_d1_pat, relevantSeqs_gr, relevantSeqs)


motif_disrupted_regions <- unique(motif_destruction_info_table_d1$region)
d1_fullGTEX_inOurSet_candidates_withMotifs <- c()
for (i in 1:length(motif_disrupted_regions)) {
  thisSet_candidates <- d1_fullGTEX_inOurSet_candidates[d1_fullGTEX_inOurSet_candidates$tf_region==motif_disrupted_regions[i],]
  this_motif_disrupted <- motif_destruction_info_table_d1[motif_destruction_info_table_d1$region==motif_disrupted_regions[i],]

  thisSet_candidates <- unique(thisSet_candidates[,c("alt_variant_id", "gene_name", "rna_region", "tf_region", "tf_names", "tf_pval", "tf_refReads", "tf_altReads", "rna_refReads", "rna_altReads", "log_rna_ratio", "gene_chr", "gene_start", "gene_end", "strand", "slope")])

  thisSet_candidates <- thisSet_candidates[thisSet_candidates$tf_names!="",]

  if(nrow(thisSet_candidates)>0) {
    allFoundMotifs <- unique(unlist(strsplit(thisSet_candidates$tf_names, split=";")))
    this_motif_disrupted_rightMotifs <- this_motif_disrupted[this_motif_disrupted$motif_alt_id%in%allFoundMotifs,]
    check <- merge(thisSet_candidates, this_motif_disrupted_rightMotifs, by.x="tf_region", by.y="region")
    d1_fullGTEX_inOurSet_candidates_withMotifs <- rbind(d1_fullGTEX_inOurSet_candidates_withMotifs, check)
  }
}

eQTL_locs_d1 <- as.data.frame(matrix(nrow=nrow(d1_fullGTEX_inOurSet_candidates_withMotifs), ncol=4, data=unlist(strsplit(d1_fullGTEX_inOurSet_candidates_withMotifs$alt_variant_id, split="_")), byrow=T))
eQTL_locs_d1 <- eQTL_locs_d1[,c(1,2,2)]
colnames(eQTL_locs_d1) <- c("chr", "start", "end")
eQTL_locs_d1_gr <- makeGRangesFromDataFrame(eQTL_locs_d1)

motif_locs_d1 <- d1_fullGTEX_inOurSet_candidates_withMotifs[,c("seqnames", "start", "end")]
motif_locs_d1_gr <- makeGRangesFromDataFrame(motif_locs_d1)

theIntersect <- as.data.frame(findOverlaps(eQTL_locs_d1_gr, motif_locs_d1_gr))
theIntersect <- theIntersect[theIntersect[,1]==theIntersect[,2],]
d1_fullGTEX_inOurSet_candidates_withMotifs <- d1_fullGTEX_inOurSet_candidates_withMotifs[unique(theIntersect[,1]),]


fimo_d2_mat <- read.table(fimo_d2_mat, header=T, sep="\t", stringsAsFactors=F)
fimo_d2_pat <- read.table(fimo_d2_pat, header=T, sep="\t", stringsAsFactors=F)

relevantSeqs <- matrix(nrow=nrow(d2_fullGTEX_inOurSet_candidates), ncol=3, data=unlist(strsplit(d2_fullGTEX_inOurSet_candidates$tf_region, split=":|-")), byrow=T)
relevantSeqs <- unique(as.data.frame(relevantSeqs))
colnames(relevantSeqs) <- c("chr", "start", "end")
relevantSeqs_gr <- makeGRangesFromDataFrame(relevantSeqs)
relevantSeqs_gr <- unique(relevantSeqs_gr)

motif_destruction_info_table_d2 <- determine_motif_destruction(fimo_d2_mat, fimo_d2_pat, relevantSeqs_gr, relevantSeqs)


motif_disrupted_regions <- unique(motif_destruction_info_table_d2$region)
d2_fullGTEX_inOurSet_candidates_withMotifs <- c()
for (i in 1:length(motif_disrupted_regions)) {
  thisSet_candidates <- d2_fullGTEX_inOurSet_candidates[d2_fullGTEX_inOurSet_candidates$tf_region==motif_disrupted_regions[i],]
  this_motif_disrupted <- motif_destruction_info_table_d2[motif_destruction_info_table_d2$region==motif_disrupted_regions[i],]

  thisSet_candidates <- unique(thisSet_candidates[,c("alt_variant_id", "gene_name", "rna_region", "tf_region", "tf_names", "tf_pval", "tf_refReads", "tf_altReads", "rna_refReads", "rna_altReads", "log_rna_ratio", "gene_chr", "gene_start", "gene_end", "strand", "slope")])

  unique_motifs <- unique(this_motif_disrupted$motif_alt_id)
  thisSet_candidates <- thisSet_candidates[thisSet_candidates$tf_names!="",]

  if(nrow(thisSet_candidates)>0) {
    allFoundMotifs <- unique(unlist(strsplit(thisSet_candidates$tf_names, split=";")))
    this_motif_disrupted_rightMotifs <- this_motif_disrupted[this_motif_disrupted$motif_alt_id%in%allFoundMotifs,]
    check <- merge(thisSet_candidates, this_motif_disrupted_rightMotifs, by.x="tf_region", by.y="region")
    d2_fullGTEX_inOurSet_candidates_withMotifs <- rbind(d2_fullGTEX_inOurSet_candidates_withMotifs, check)
  }
}


eQTL_locs_d2 <- as.data.frame(matrix(nrow=nrow(d2_fullGTEX_inOurSet_candidates_withMotifs), ncol=4, data=unlist(strsplit(d2_fullGTEX_inOurSet_candidates_withMotifs$alt_variant_id, split="_")), byrow=T))
eQTL_locs_d2 <- eQTL_locs_d2[,c(1,2,2)]
colnames(eQTL_locs_d2) <- c("chr", "start", "end")
eQTL_locs_d2_gr <- makeGRangesFromDataFrame(eQTL_locs_d2)

motif_locs_d2 <- d2_fullGTEX_inOurSet_candidates_withMotifs[,c("seqnames", "start", "end")]
motif_locs_d2_gr <- makeGRangesFromDataFrame(motif_locs_d2)

theIntersect <- as.data.frame(findOverlaps(eQTL_locs_d2_gr, motif_locs_d2_gr))
theIntersect <- theIntersect[theIntersect[,1]==theIntersect[,2],]
d2_fullGTEX_inOurSet_candidates_withMotifs <- d2_fullGTEX_inOurSet_candidates_withMotifs[unique(theIntersect[,1]),]





################################################################################
#Define outdir to which plots will be written.
################################################################################

plot_outDir <- paste(outDir, "plotgardener_promoter_motif/", sep="")
if(!file.exists(plot_outDir)) {system(paste("mkdir ", plot_outDir))}


################################################################################
#Read in HiC and make appropriate genomic ranges.
################################################################################


HiC <- "/cluster/home/lrizzardi/HiC_analysis/hiccups_220622/GLUTA_HiC_5kb_full/merged_loops.bedpe"
HiC <- read.table(HiC, header=F, sep="\t", stringsAsFactors=F)

HiC_1_gr <- as.data.frame(HiC[,1:3])
HiC_1_gr[,1] <- paste("chr", HiC_1_gr[,1], sep="")
colnames(HiC_1_gr) <- c("chr", "start", "end")
HiC_1_gr <- makeGRangesFromDataFrame(HiC_1_gr)

HiC_2_gr <- as.data.frame(HiC[,4:6])
HiC_2_gr[,1] <- paste("chr", HiC_2_gr[,1], sep="")
colnames(HiC_2_gr) <- c("chr", "start", "end")
HiC_2_gr <- makeGRangesFromDataFrame(HiC_2_gr)



################################################################################
#For each of these, use plotgardener to just plot
#1) All the relevant gene regions,
#2) The noted connections in the HiC data
#3) The connections in our dataset.
#4) Make barplots for the RNA and TF reference and alt reads.
################################################################################



################################################################################
#Plots for donor 1.
################################################################################

unique_regions <- unique(d1_fullGTEX_inOurSet_candidates_withMotifs[,1])

for (i in 1:length(unique_regions)) {
  print(paste(i, length(unique_regions)))
  saveBase <- paste(plot_outDir, "d1_", gsub("-", "_", gsub(":", "_", unique_regions[i])), sep="")

  thisSet <- d1_fullGTEX_inOurSet_candidates_withMotifs[d1_fullGTEX_inOurSet_candidates_withMotifs[,1]==unique_regions[i],]

  theSlopes <- unique(thisSet$slope)
  theSlopes <- paste(theSlopes, collapse=", ")
  eQTL_locs <- as.numeric(unlist(strsplit(thisSet[,2], split="_")))
  eQTL_locs <- unique(eQTL_locs[!is.na(eQTL_locs)])
  all_motifs <- paste(unique(paste(thisSet$motif_id, thisSet$strand.y, thisSet$start, thisSet$end, sep="_")), collapse=", ")
  all_motifs_alt <- paste(unique(thisSet$motif_alt_id), collapse=", ")
  eQTL_label <- paste(unique(thisSet[,2]), collapse=", ")

  thisSet <- unique(thisSet[,c("alt_variant_id", "gene_chr", "gene_start", "gene_end", "strand.x", "tf_region", "rna_region", "rna_refReads", "rna_altReads", "tf_refReads", "tf_altReads", "tf_names", "gene_name")])


  ###
  #Identify all TFs affected, their p-values, etc
  ###
  all_TFs <- unlist(strsplit(thisSet[1,"tf_names"], split=";"))
  all_TFs_refReads <- as.numeric(unlist(strsplit(thisSet[1,"tf_refReads"], split=";")))
  all_TFs_altReads <- as.numeric(unlist(strsplit(thisSet[1,"tf_altReads"], split=";")))
  all_TFs_logRatio <- log((all_TFs_refReads+1)/(all_TFs_altReads+1))

  ###
  #Make a plot of all of the TF reads.  First, do raw read counts.
  ###
  graphDF_1 <- data.frame(TF=all_TFs, Reads=all_TFs_refReads, type=rep("Reference", length(all_TFs)))
  graphDF_2 <- data.frame(TF=all_TFs, Reads=all_TFs_altReads, type=rep("Alternate", length(all_TFs)))
  graphDF <- rbind(graphDF_1, graphDF_2)
  graphDF$TF <- factor(graphDF$TF, levels=all_TFs[order(all_TFs)])

  saveFile_1 <- paste(saveBase, "allTFs_reads_barplot.pdf", sep="")
  p_tf_barplot <- ggplot(data=graphDF, aes(x=TF, y=Reads, fill=type)) +
    geom_bar(stat="identity", position="dodge") + theme_classic(base_size=18) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="bottom") + ylab("ChIP-seq Reads")
  ggsave(saveFile_1)


  ###
  #Make a plot of all of the RNA reads.  First, do raw read counts.
  ###
  graphDF_1 <- unique(data.frame(Regions=thisSet$rna_region, Reads=thisSet$rna_refReads, type=rep("Reference", nrow(thisSet))))
  graphDF_1$Region <- c(1:nrow(graphDF_1))
  graphDF_2 <- unique(data.frame(Regions=thisSet$rna_region, Reads=thisSet$rna_altReads, type=rep("Alternate", nrow(thisSet))))
  graphDF_2$Region <- c(1:nrow(graphDF_2))
  graphDF <- rbind(graphDF_1, graphDF_2)
  graphDF <- unique(graphDF)
  rna_regions <- unique(graphDF$Region)
  rna_regions <- rna_regions[order(rna_regions)]
  graphDF$Region <- factor(graphDF$Region, levels=rna_regions)

  saveFile_1 <- paste(saveBase, "allRNA_reads_barplot.pdf", sep="")
  p_tf_barplot <- ggplot(data=graphDF, aes(x=Region, y=Reads, fill=type)) +
    geom_bar(stat="identity", position="dodge") + theme_classic(base_size=18) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="bottom") + ylab("RNA Reads")
  ggsave(saveFile_1)


  tf_regions <- unique(as.data.frame(matrix(nrow=nrow(thisSet), ncol=3, data=unlist(strsplit(thisSet$tf_region, split=":|-")), byrow=T)))
  tf_regions[,2] <- as.numeric(tf_regions[,2])
  tf_regions[,3] <- as.numeric(tf_regions[,3])

  rna_regions <- unique(as.data.frame(matrix(nrow=nrow(thisSet), ncol=3, data=unlist(strsplit(thisSet$rna_region, split=":|-")), byrow=T)))
  rna_regions[,2] <- as.numeric(rna_regions[,2])
  rna_regions[,3] <- as.numeric(rna_regions[,3])
  rna_regions$locations <- rowMeans(as.matrix(rna_regions[,2:3]))
  rna_regions$label <- c(1:nrow(rna_regions))

  theChr <- thisSet[1,"gene_chr"]
  theLocs <- c(unique(tf_regions[,2]), unique(tf_regions[,3]), unique(rna_regions[,2]), unique(rna_regions[,3]), unique(thisSet[,"gene_start"]), unique(thisSet[,"gene_end"]))
  theLocs <- as.numeric(theLocs)
  chrMin <- min(theLocs)-2000
  chrMax <- max(theLocs)+2000

  new_gr <- makeGRangesFromDataFrame(data.frame(chr=theChr, start=min(eQTL_locs), end=max(eQTL_locs)))
  interSect1 <- as.data.frame(findOverlaps(new_gr, HiC_1_gr))
  interSect2 <- as.data.frame(findOverlaps(new_gr, HiC_2_gr))
  allHiCOverlaps <- unique(c(interSect1[,2], interSect2[,2]))

  saveFile <- paste(saveBase, ".pdf", sep="")
  pdf(saveFile)
  pageCreate(width = 8, height = 8.5, default.units = "inches")
  params_a <- pgParams(chrom = theChr, chromstart = chrMin, chromend = chrMax, assembly = "hg38",
                x = 0.75, width = 5, just = c("left", "top"), default.units = "inches")
  #
  if(length(allHiCOverlaps)>0) {
    this_HiC_PE <- HiC[allHiCOverlaps,]
    alt_theLocs <- unique(as.numeric(c(this_HiC_PE[,2], this_HiC_PE[,3], this_HiC_PE[,5], this_HiC_PE[,6])))
    alt_chrMin <- min(c(theLocs, alt_theLocs))-2000
    alt_chrMax <- max(c(theLocs, alt_theLocs))+2000
    this_HiC_PE <- this_HiC_PE[,1:6]
    this_HiC_PE[,1] <- paste("chr", this_HiC_PE[,1], sep="")
    this_HiC_PE[,4] <- paste("chr", this_HiC_PE[,4], sep="")
    colnames(this_HiC_PE) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    this_HiC_PE$height <- rep(1, nrow(this_HiC_PE))

    params_b <- pgParams(chrom = theChr, chromstart = alt_chrMin, chromend = alt_chrMax, assembly = "hg38",
                  x = 0.75, width = 5, just = c("left", "top"), default.units = "inches")
    #
    bedpePlot_d2 <- plotPairsArches(data = this_HiC_PE, params = params_b,
        archHeight = this_HiC_PE$height, alpha = 1, color="black", height = 0.5, just = c("left", "top"), y=1.25,
        default.units = "inches")
    genes_b <- plotGenes(params = params_b, stroke = 1, fontsize = 10, y = 2, height = 0.4)
    annoGenomeLabel(plot = genes_b, params = params_a, scale = "Kb", fontsize = 7, y = 1.75)
  }


  genes_a <- plotGenes(params = params_a, stroke = 1, fontsize = 10, y = 4, height = 0.4)
  annoGenomeLabel(plot = genes_a, params = params_a, scale = "Kb", fontsize = 7, y = 3.75)
  annoText(label="*", plot=genes_a, x=eQTL_locs, y= -0.15, default.units="native", fontsize = 10)
  annoText(label="|", plot=genes_a, x=rna_regions$locations, y= -0.15, default.units="native", fontsize = 10)
  annoText(label=rna_regions$label, plot=genes_a, x=rna_regions$locations, y= -0.4, default.units="native", fontsize = 10)
  annoText(label=all_motifs, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -1.5, default.units="native", fontsize = 5)
  annoText(label=theSlopes, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -2, default.units="native", fontsize = 5)
  annoText(label=all_motifs_alt, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -2.5, default.units="native", fontsize = 5)
  annoText(label=eQTL_label, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -3, default.units="native", fontsize = 5)
  #annoText(label=all_motifs_alt, plot=genes_a, x=eQTL_locs[1], y= 3.5, default.units="native", fontsize = 10)
  #annoText(label=eQTL_label, plot=genes_a, x=eQTL_locs[1], y= 4.5, default.units="native", fontsize = 10)

  pageGuideHide()
  dev.off()

}




################################################################################
#Plots for donor 2.
################################################################################

unique_regions <- unique(d2_fullGTEX_inOurSet_candidates_withMotifs[,1])

for (i in 1:length(unique_regions)) {
  print(paste(i, length(unique_regions)))
  saveBase <- paste(plot_outDir, "d2_", gsub("-", "_", gsub(":", "_", unique_regions[i])), sep="")

  thisSet <- d2_fullGTEX_inOurSet_candidates_withMotifs[d2_fullGTEX_inOurSet_candidates_withMotifs[,1]==unique_regions[i],]

  theSlopes <- unique(thisSet$slope)
  theSlopes <- paste(theSlopes, collapse=", ")
  eQTL_locs <- as.numeric(unlist(strsplit(thisSet[,2], split="_")))
  eQTL_locs <- unique(eQTL_locs[!is.na(eQTL_locs)])
  all_motifs <- paste(unique(paste(thisSet$motif_id, thisSet$strand.y, thisSet$start, thisSet$end, sep="_")), collapse=", ")
  all_motifs_alt <- paste(unique(thisSet$motif_alt_id), collapse=", ")
  eQTL_label <- paste(unique(thisSet[,2]), collapse=", ")

  thisSet <- unique(thisSet[,c("alt_variant_id", "gene_chr", "gene_start", "gene_end", "strand.x", "tf_region", "rna_region", "rna_refReads", "rna_altReads", "tf_refReads", "tf_altReads", "tf_names", "gene_name")])


  ###
  #Identify all TFs affected, their p-values, etc
  ###
  all_TFs <- unlist(strsplit(thisSet[1,"tf_names"], split=";"))
  all_TFs_refReads <- as.numeric(unlist(strsplit(thisSet[1,"tf_refReads"], split=";")))
  all_TFs_altReads <- as.numeric(unlist(strsplit(thisSet[1,"tf_altReads"], split=";")))
  all_TFs_logRatio <- log((all_TFs_refReads+1)/(all_TFs_altReads+1))

  ###
  #Make a plot of all of the TF reads.  First, do raw read counts.
  ###
  graphDF_1 <- data.frame(TF=all_TFs, Reads=all_TFs_refReads, type=rep("Reference", length(all_TFs)))
  graphDF_2 <- data.frame(TF=all_TFs, Reads=all_TFs_altReads, type=rep("Alternate", length(all_TFs)))
  graphDF <- rbind(graphDF_1, graphDF_2)
  graphDF$TF <- factor(graphDF$TF, levels=all_TFs[order(all_TFs)])

  saveFile_1 <- paste(saveBase, "allTFs_reads_barplot.pdf", sep="")
  p_tf_barplot <- ggplot(data=graphDF, aes(x=TF, y=Reads, fill=type)) +
    geom_bar(stat="identity", position="dodge") + theme_classic(base_size=18) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="bottom") + ylab("ChIP-seq Reads")
  ggsave(saveFile_1)


  ###
  #Make a plot of all of the RNA reads.  First, do raw read counts.
  ###
  graphDF_1 <- unique(data.frame(Regions=thisSet$rna_region, Reads=thisSet$rna_refReads, type=rep("Reference", nrow(thisSet))))
  graphDF_1$Region <- c(1:nrow(graphDF_1))
  graphDF_2 <- unique(data.frame(Regions=thisSet$rna_region, Reads=thisSet$rna_altReads, type=rep("Alternate", nrow(thisSet))))
  graphDF_2$Region <- c(1:nrow(graphDF_2))
  graphDF <- rbind(graphDF_1, graphDF_2)
  graphDF <- unique(graphDF)
  rna_regions <- unique(graphDF$Region)
  rna_regions <- rna_regions[order(rna_regions)]
  graphDF$Region <- factor(graphDF$Region, levels=rna_regions)

  saveFile_1 <- paste(saveBase, "allRNA_reads_barplot.pdf", sep="")
  p_tf_barplot <- ggplot(data=graphDF, aes(x=Region, y=Reads, fill=type)) +
    geom_bar(stat="identity", position="dodge") + theme_classic(base_size=18) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position="bottom") + ylab("RNA Reads")
  ggsave(saveFile_1)


  tf_regions <- unique(as.data.frame(matrix(nrow=nrow(thisSet), ncol=3, data=unlist(strsplit(thisSet$tf_region, split=":|-")), byrow=T)))
  tf_regions[,2] <- as.numeric(tf_regions[,2])
  tf_regions[,3] <- as.numeric(tf_regions[,3])

  rna_regions <- unique(as.data.frame(matrix(nrow=nrow(thisSet), ncol=3, data=unlist(strsplit(thisSet$rna_region, split=":|-")), byrow=T)))
  rna_regions[,2] <- as.numeric(rna_regions[,2])
  rna_regions[,3] <- as.numeric(rna_regions[,3])
  rna_regions$locations <- rowMeans(as.matrix(rna_regions[,2:3]))
  rna_regions$label <- c(1:nrow(rna_regions))

  theChr <- thisSet[1,"gene_chr"]
  theLocs <- c(unique(tf_regions[,2]), unique(tf_regions[,3]), unique(rna_regions[,2]), unique(rna_regions[,3]), unique(thisSet[,"gene_start"]), unique(thisSet[,"gene_end"]))
  theLocs <- as.numeric(theLocs)
  chrMin <- min(theLocs)-2000
  chrMax <- max(theLocs)+2000

  new_gr <- makeGRangesFromDataFrame(data.frame(chr=theChr, start=min(eQTL_locs), end=max(eQTL_locs)))
  interSect1 <- as.data.frame(findOverlaps(new_gr, HiC_1_gr))
  interSect2 <- as.data.frame(findOverlaps(new_gr, HiC_2_gr))
  allHiCOverlaps <- unique(c(interSect1[,2], interSect2[,2]))

  saveFile <- paste(saveBase, ".pdf", sep="")
  pdf(saveFile)
  pageCreate(width = 8, height = 8.5, default.units = "inches")
  params_a <- pgParams(chrom = theChr, chromstart = chrMin, chromend = chrMax, assembly = "hg38",
                x = 0.75, width = 5, just = c("left", "top"), default.units = "inches")
  #
  if(length(allHiCOverlaps)>0) {
    this_HiC_PE <- HiC[allHiCOverlaps,]
    alt_theLocs <- unique(as.numeric(c(this_HiC_PE[,2], this_HiC_PE[,3], this_HiC_PE[,5], this_HiC_PE[,6])))
    alt_chrMin <- min(c(theLocs, alt_theLocs))-2000
    alt_chrMax <- max(c(theLocs, alt_theLocs))+2000
    this_HiC_PE <- this_HiC_PE[,1:6]
    this_HiC_PE[,1] <- paste("chr", this_HiC_PE[,1], sep="")
    this_HiC_PE[,4] <- paste("chr", this_HiC_PE[,4], sep="")
    colnames(this_HiC_PE) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    this_HiC_PE$height <- rep(1, nrow(this_HiC_PE))

    params_b <- pgParams(chrom = theChr, chromstart = alt_chrMin, chromend = alt_chrMax, assembly = "hg38",
                  x = 0.75, width = 5, just = c("left", "top"), default.units = "inches")
    #
    bedpePlot_d2 <- plotPairsArches(data = this_HiC_PE, params = params_b,
        archHeight = this_HiC_PE$height, alpha = 1, color="black", height = 0.5, just = c("left", "top"), y=1.25,
        default.units = "inches")
    genes_b <- plotGenes(params = params_b, stroke = 1, fontsize = 10, y = 2, height = 0.4)
    annoGenomeLabel(plot = genes_b, params = params_a, scale = "Kb", fontsize = 7, y = 1.75)
  }


  genes_a <- plotGenes(params = params_a, stroke = 1, fontsize = 10, y = 4, height = 0.4)
  annoGenomeLabel(plot = genes_a, params = params_a, scale = "Kb", fontsize = 7, y = 3.75)
  annoText(label="*", plot=genes_a, x=eQTL_locs, y= -0.15, default.units="native", fontsize = 10)
  annoText(label="|", plot=genes_a, x=rna_regions$locations, y= -0.15, default.units="native", fontsize = 10)
  annoText(label=rna_regions$label, plot=genes_a, x=rna_regions$locations, y= -0.4, default.units="native", fontsize = 10)
  annoText(label=all_motifs, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -1.5, default.units="native", fontsize = 5)
  annoText(label=theSlopes, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -2, default.units="native", fontsize = 5)
  annoText(label=all_motifs_alt, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -2.5, default.units="native", fontsize = 5)
  annoText(label=eQTL_label, plot=genes_a, x=mean(c(chrMin, chrMax)), y= -3, default.units="native", fontsize = 5)
  #annoText(label=all_motifs_alt, plot=genes_a, x=eQTL_locs[1], y= 3.5, default.units="native", fontsize = 10)
  #annoText(label=eQTL_label, plot=genes_a, x=eQTL_locs[1], y= 4.5, default.units="native", fontsize = 10)

  pageGuideHide()
  dev.off()

}






################################################################################
################################################################################
#We now want to make a few bespoke plots for specific examples, based on the
#above results.  In particular, we want to make plots which show the
#FRACTION of reads which map to the reference and alternate alleles in each case
#for the FBXO7 example (shared in both donors) and the RPS14 example
#(found only in donor 2).
#For the RNA case, we want to sum reads across all relevant het regions,
#rather than showing them separately.
#For the TF case, we will of course display each TF separately.
#Finally, for the FBX07 case, we also want to determine the reads mapping to
#SYN3 as well, since the eQTL of interest also happens to be
#an eQTL for SYN3, with a negative slope rather than a positive one.
################################################################################
################################################################################


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

geneBodies <- GenomicFeatures::genes(txdb)
#8224 is SYN3
#100130890 is TSTD3
#85015 is USP45
SYN3_gr <- geneBodies[(elementMetadata(geneBodies)[, "gene_id"]==8224)]
FBX07_gr <- geneBodies[(elementMetadata(geneBodies)[, "gene_id"]==25793)]
RPS14_gr <- geneBodies[(elementMetadata(geneBodies)[, "gene_id"]==6208)]
TSTD3_gr <- geneBodies[(elementMetadata(geneBodies)[, "gene_id"]==100130890)]
USP45_gr <- geneBodies[(elementMetadata(geneBodies)[, "gene_id"]==85015)]

altDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/"
saveFile <- paste(altDir, "rna_sig_d1_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
d1_RNA_TF_linked <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)


saveFile <- paste(altDir, "rna_sig_d2_regions_wPhase_w_SignificantTFInPhase.txt", sep="")
d2_RNA_TF_linked <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)




################################################################################
#FBX07
#In this case, maternal is the alternate, paternal is the reference
################################################################################

d1_RNA_TF_linked$tf_region <- paste(d1_RNA_TF_linked[,"chr_tf"], ":", d1_RNA_TF_linked[,"start_tf"], "-", d1_RNA_TF_linked[,"end_tf"], sep="")
relevant_d1 <- d1_RNA_TF_linked[d1_RNA_TF_linked$tf_region=="chr22:32474682-32474882",]
relevant_d1_gr <- makeGRangesFromDataFrame(relevant_d1)
theIntersect <- as.data.frame(findOverlaps(relevant_d1_gr, SYN3_gr))
relevant_d1_SYN3 <- relevant_d1[unique(theIntersect[,1]),]
theIntersect <- as.data.frame(findOverlaps(relevant_d1_gr, FBX07_gr))
relevant_d1_FBX07 <- relevant_d1[unique(theIntersect[,1]),]



d2_RNA_TF_linked$tf_region <- paste(d2_RNA_TF_linked[,"chr_tf"], ":", d2_RNA_TF_linked[,"start_tf"], "-", d2_RNA_TF_linked[,"end_tf"], sep="")
relevant_d2 <- d2_RNA_TF_linked[d2_RNA_TF_linked$tf_region=="chr22:32474682-32474919",]
relevant_d2_gr <- makeGRangesFromDataFrame(relevant_d2)
theIntersect <- as.data.frame(findOverlaps(relevant_d2_gr, SYN3_gr))
relevant_d2_SYN3 <- relevant_d2[unique(theIntersect[,1]),]
theIntersect <- as.data.frame(findOverlaps(relevant_d2_gr, FBX07_gr))
relevant_d2_FBX07 <- relevant_d2[unique(theIntersect[,1]),]


################################################################################
#From inspection of the above graphs, I know that Mat is alternate and
#pat is reference for this phased block.
#Set up the graphDF.
################################################################################

summed_RNA_reads_FBX07_d1 <- as.data.frame(rbind(colSums(relevant_d1_FBX07[,c("rna_matReads", "rna_patReads")])))
summed_RNA_reads_FBX07_d1$Donor <- rep("Donor1", nrow(summed_RNA_reads_FBX07_d1))
summed_RNA_reads_FBX07_d2 <- as.data.frame(rbind(colSums(relevant_d2_FBX07[,c("rna_matReads", "rna_patReads")])))
summed_RNA_reads_FBX07_d2$Donor <- rep("Donor2", nrow(summed_RNA_reads_FBX07_d2))
summed_RNA_reads_FBX07 <- rbind(summed_RNA_reads_FBX07_d1, summed_RNA_reads_FBX07_d2)
summed_RNA_reads_FBX07$gene <- rep("FBXO7", nrow(summed_RNA_reads_FBX07))
colnames(summed_RNA_reads_FBX07)[1:2] <- c("Alternate", "Reference")
summed_RNA_reads_FBX07$AlternateFraction <- summed_RNA_reads_FBX07$Alternate / rowSums(summed_RNA_reads_FBX07[,c("Alternate", "Reference")])
summed_RNA_reads_FBX07$ReferenceFraction <- summed_RNA_reads_FBX07$Reference / rowSums(summed_RNA_reads_FBX07[,c("Alternate", "Reference")])

summed_RNA_reads_SYN3_d1 <- as.data.frame(rbind(colSums(relevant_d1_SYN3[,c("rna_matReads", "rna_patReads")])))
summed_RNA_reads_SYN3_d1$Donor <- rep("Donor1", nrow(summed_RNA_reads_SYN3_d1))
summed_RNA_reads_SYN3_d2 <- as.data.frame(rbind(colSums(relevant_d2_SYN3[,c("rna_matReads", "rna_patReads")])))
summed_RNA_reads_SYN3_d2$Donor <- rep("Donor2", nrow(summed_RNA_reads_SYN3_d2))
summed_RNA_reads_SYN3 <- rbind(summed_RNA_reads_SYN3_d1, summed_RNA_reads_SYN3_d2)
summed_RNA_reads_SYN3$gene <- rep("SYN3", nrow(summed_RNA_reads_SYN3))
colnames(summed_RNA_reads_SYN3)[1:2] <- c("Alternate", "Reference")
summed_RNA_reads_SYN3$AlternateFraction <- summed_RNA_reads_SYN3$Alternate / rowSums(summed_RNA_reads_SYN3[,c("Alternate", "Reference")])
summed_RNA_reads_SYN3$ReferenceFraction <- summed_RNA_reads_SYN3$Reference / rowSums(summed_RNA_reads_SYN3[,c("Alternate", "Reference")])


graphDF <- rbind(summed_RNA_reads_FBX07, summed_RNA_reads_SYN3)
graphDF_1 <- graphDF[,c(1:5)]
colnames(graphDF_1)[5] <- "Fraction"
graphDF_1$Haplotype <- rep("Alternate", nrow(graphDF_1))
graphDF_2 <- graphDF[,c(1:4,6)]
colnames(graphDF_2)[5] <- "Fraction"
graphDF_2$Haplotype <- rep("Reference", nrow(graphDF_1))
graphDF <- rbind(graphDF_1, graphDF_2)

saveFile <- paste(plot_outDir, "Example_FBXO7_SYN3_RNA_fraction_Stacked_Barplot.pdf", sep="")
p_tf_barplot <- ggplot(data=graphDF, aes(x=gene, y=Fraction, fill=Haplotype)) + facet_wrap(~Donor) +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="bottom") + ylab("RNA Reads Fraction")
ggsave(saveFile)




################################################################################
#The supplemental will have information for all factors.  But let's go ahead
#and make the main figure have information for just CTCF and RAD21, specifically.
################################################################################


donor1_TFs <- data.frame(TF=unlist(strsplit(relevant_d1_FBX07[1,"tf_names"], split=";")),
  Alternate=unlist(strsplit(relevant_d1_FBX07[1,"tf_matReads"], split=";")),
  Reference=unlist(strsplit(relevant_d1_FBX07[1,"tf_patReads"], split=";")))
donor1_TFs$Alternate <- as.numeric(as.character(donor1_TFs$Alternate))
donor1_TFs$Reference <- as.numeric(as.character(donor1_TFs$Reference))

donor2_TFs <- data.frame(TF=unlist(strsplit(relevant_d2_FBX07[1,"tf_names"], split=";")),
  Alternate=unlist(strsplit(relevant_d2_FBX07[1,"tf_matReads"], split=";")),
  Reference=unlist(strsplit(relevant_d2_FBX07[1,"tf_patReads"], split=";")))
donor2_TFs$Alternate <- as.numeric(as.character(donor2_TFs$Alternate))
donor2_TFs$Reference <- as.numeric(as.character(donor2_TFs$Reference))

#donor1_TFs <- donor1_TFs[donor1_TFs$TF%in%donor2_TFs$TF,]
#donor2_TFs <- donor2_TFs[donor2_TFs$TF%in%donor1_TFs$TF,]

donor1_TFs$Donor <- rep("Donor1", nrow(donor1_TFs))
donor2_TFs$Donor <- rep("Donor2", nrow(donor2_TFs))

donor1_TFs$AlternateFraction <- donor1_TFs$Alternate / rowSums(donor1_TFs[,c("Alternate", "Reference")])
donor1_TFs$ReferenceFraction <- donor1_TFs$Reference / rowSums(donor1_TFs[,c("Alternate", "Reference")])

donor2_TFs$AlternateFraction <- donor2_TFs$Alternate / rowSums(donor2_TFs[,c("Alternate", "Reference")])
donor2_TFs$ReferenceFraction <- donor2_TFs$Reference / rowSums(donor2_TFs[,c("Alternate", "Reference")])


graphDF <- rbind(donor1_TFs, donor2_TFs)
graphDF_1 <- graphDF[,c(1:5)]
colnames(graphDF_1)[5] <- "Fraction"
graphDF_1$Haplotype <- rep("Alternate", nrow(graphDF_1))
graphDF_2 <- graphDF[,c(1:4,6)]
colnames(graphDF_2)[5] <- "Fraction"
graphDF_2$Haplotype <- rep("Reference", nrow(graphDF_1))
graphDF <- rbind(graphDF_1, graphDF_2)


saveFile <- paste(plot_outDir, "Example_FBXO7_SYN3_TF_fraction_Stacked_Barplot.pdf", sep="")
p_tf_barplot <- ggplot(data=graphDF, aes(x=TF, y=Fraction, fill=Haplotype)) + facet_wrap(~Donor) +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="bottom") + ylab("ChIP-seq Reads Fraction")
ggsave(saveFile)




################################################################################
#RSP14
#In this case, maternal is the reference, paternal is the alternate.
################################################################################

relevant_d2 <- d2_RNA_TF_linked[d2_RNA_TF_linked$tf_region=="chr5:150449648-150449848",]
relevant_d2_gr <- makeGRangesFromDataFrame(relevant_d2)
theIntersect <- as.data.frame(findOverlaps(relevant_d2_gr, RPS14_gr))
relevant_d2_RPS14 <- relevant_d2[unique(theIntersect[,1]),]


################################################################################
#From inspection of the above graphs, I know that Mat is reference and
#pat is alternate for this phased block.
#Set up the graphDF.
################################################################################

summed_RNA_reads_RPS14 <- as.data.frame(rbind(colSums(relevant_d2_RPS14[,c("rna_matReads", "rna_patReads")])))
summed_RNA_reads_RPS14$Donor <- rep("Donor2", nrow(summed_RNA_reads_RPS14))
summed_RNA_reads_RPS14$gene <- rep("RPS14", nrow(summed_RNA_reads_RPS14))
colnames(summed_RNA_reads_RPS14)[1:2] <- c("Reference", "Alternate")
summed_RNA_reads_RPS14$AlternateFraction <- summed_RNA_reads_RPS14$Alternate / rowSums(summed_RNA_reads_RPS14[,c("Alternate", "Reference")])
summed_RNA_reads_RPS14$ReferenceFraction <- summed_RNA_reads_RPS14$Reference / rowSums(summed_RNA_reads_RPS14[,c("Alternate", "Reference")])


graphDF <- rbind(summed_RNA_reads_RPS14)
graphDF_1 <- graphDF[,c(1:5)]
colnames(graphDF_1)[5] <- "Fraction"
graphDF_1$Haplotype <- rep("Alternate", nrow(graphDF_1))
graphDF_2 <- graphDF[,c(1:4,6)]
colnames(graphDF_2)[5] <- "Fraction"
graphDF_2$Haplotype <- rep("Reference", nrow(graphDF_1))
graphDF <- rbind(graphDF_1, graphDF_2)

saveFile <- paste(plot_outDir, "Example_RPS14_RNA_fraction_Stacked_Barplot.pdf", sep="")
p_tf_barplot <- ggplot(data=graphDF, aes(x=gene, y=Fraction, fill=Haplotype)) + facet_wrap(~Donor) +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="bottom") + ylab("RNA Reads Fraction")
ggsave(saveFile)




################################################################################
#Show all TFs.
################################################################################

donor2_TFs <- data.frame(TF=unlist(strsplit(relevant_d2_RPS14[1,"tf_names"], split=";")),
  Reference=unlist(strsplit(relevant_d2_RPS14[1,"tf_matReads"], split=";")),
  Alternate=unlist(strsplit(relevant_d2_RPS14[1,"tf_patReads"], split=";")))
donor2_TFs$Alternate <- as.numeric(as.character(donor2_TFs$Alternate))
donor2_TFs$Reference <- as.numeric(as.character(donor2_TFs$Reference))

donor2_TFs$Donor <- rep("Donor2", nrow(donor2_TFs))

donor2_TFs$AlternateFraction <- donor2_TFs$Alternate / rowSums(donor2_TFs[,c("Alternate", "Reference")])
donor2_TFs$ReferenceFraction <- donor2_TFs$Reference / rowSums(donor2_TFs[,c("Alternate", "Reference")])


graphDF <- rbind(donor2_TFs)
graphDF_1 <- graphDF[,c(1:5)]
colnames(graphDF_1)[5] <- "Fraction"
graphDF_1$Haplotype <- rep("Alternate", nrow(graphDF_1))
graphDF_2 <- graphDF[,c(1:4,6)]
colnames(graphDF_2)[5] <- "Fraction"
graphDF_2$Haplotype <- rep("Reference", nrow(graphDF_1))
graphDF <- rbind(graphDF_1, graphDF_2)


saveFile <- paste(plot_outDir, "Example_RPS14_TF_fraction_Stacked_Barplot.pdf", sep="")
p_tf_barplot <- ggplot(data=graphDF, aes(x=TF, y=Fraction, fill=Haplotype)) + facet_wrap(~Donor) +
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="bottom") + ylab("ChIP-seq Reads Fraction")
ggsave(saveFile)





#
