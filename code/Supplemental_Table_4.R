#!/usr/bin/env R
#Supplemental_Table_4.R


################################################################################
#The following script is used to produce Figure 3D and Supplemental Figure 8. It
#     looks for disruption of motifs in heterozygous sequences, and determines,
#     for each disrupted motif, the fraction of reads favoring the ancestral
#     versus derived sequence.
#
#Run under R 4.1.0, but can be used under any other version that has the following
#     packages available:
#     matrixStats
#     ggplot2
#     GenomicRanges
#     gridExtra
#     memes
#     universalmotif
#
#This program takes the following arguments:
# outDir: The directory to which results should be written.
# theMotifs: The JASPAR motif database.
# pval_d1: A table compiling p-value calculations when reads are summed across
#     experiments for donor1, produced by the script CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# pval_d2: As pval_d1, but for donor 2.
# maternal_fasta_d1: A fasta file containing all of the sequences for haplotype 1
#     in donor 1, a subset of sequences produced by the script Construct_fasta_variant_sequences.R
# paternal_fasta_d1: As maternal_fasta_d1, but for haplotype 2.
# maternal_fasta_d2: As maternal_fasta_d1, but for donor 2.
# paternal_fasta_d2: As maternal_fasta_d2, but for haplotype 2.
# suppTable1: Supplemental Table 1.
# suppTable2: Supplemental Table 2.
#
################################################################################

#module load cluster/R/4.1.0
#cd /cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/
#module load cluster/util/gcc/8.2


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

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################

################################################################################
#The following function determines density of points in a region for coloring
#in a dot plot.
################################################################################
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


################################################################################
#The following function determines changes in motif strength between alleles.
#1) First, for a given variant sequence, identify all of the motifs which were found
#   in said sequence in both the maternal and paternal.
#2) Second, identify those which were entirely removed or gained. Save in a table.
#3) Third, identify cases of the same motif which has had a change in score, and save
#   in a table.
################################################################################
determine_motif_destruction <- function(d1_maternal_Hits, d1_paternal_Hits, relevantSeqs_gr) {
  #d1_maternal_Hits <- maternal_Hits
  #d1_paternal_Hits <- paternal_Hits

  theIntersect_maternal <- as.data.frame(findOverlaps(d1_maternal_Hits, relevantSeqs_gr))
  theIntersect_paternal <- as.data.frame(findOverlaps(d1_paternal_Hits, relevantSeqs_gr))

  uniqueRegionIndices <- unique(c(unique(theIntersect_maternal[,2]), unique(theIntersect_paternal[,2])))
  uniqueRegionIndices <- uniqueRegionIndices[order(uniqueRegionIndices)]

  motif_destruction_info_table <- c()
  for (i in 1:length(uniqueRegionIndices)) {
    thisSet_maternal <- as.data.frame(d1_maternal_Hits[theIntersect_maternal[theIntersect_maternal[,2]==uniqueRegionIndices[i],1],])
    thisSet_paternal <- as.data.frame(d1_paternal_Hits[theIntersect_paternal[theIntersect_paternal[,2]==uniqueRegionIndices[i],1],])

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
            thisLine$region <- relevantSeqs[lineIntersect[1,2]]

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
            thisLine$region <- relevantSeqs[lineIntersect[1,2]]

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
            thisLine$region <- relevantSeqs[lineIntersect[1,2]]

            motif_destruction_info_table <- rbind(motif_destruction_info_table, thisLine)
          }
        }
      }

      if(nrow(thisSet_maternal)==0 && nrow(thisSet_paternal)>0) {
        ##############################################################################
        #This represents the cases where there was a destruction of the motif
        #in maternal, but paternal has the motif.
        ##############################################################################
        for (j in 1:nrow(thisSet_paternal)) {
          thisLine <- thisSet_paternal[j,1:7]
          thisLine$score_maternal <- NA
          thisLine$pvalue_maternal <- NA
          thisLine$score_paternal <- thisSet_paternal[j,8]
          thisLine$pvalue_paternal <- thisSet_paternal[j,9]

          thisLine_gr <- makeGRangesFromDataFrame(thisLine)
          lineIntersect <- as.data.frame(findOverlaps(thisLine_gr,relevantSeqs_gr))
          thisLine$region <- relevantSeqs[lineIntersect[1,2]]

          motif_destruction_info_table <- rbind(motif_destruction_info_table, thisLine)
        }
      }
      if(nrow(thisSet_maternal)>0 && nrow(thisSet_paternal)==0) {
        ##############################################################################
        #This represents the cases where there was a destruction of the motif
        #in paternal, but maternal has the motif.
        ##############################################################################
        if(nrow(thisSet_maternal)>0) {
          for (j in 1:nrow(thisSet_maternal)) {
            thisLine <- thisSet_maternal[j,1:9]
            colnames(thisLine)[8:9] <- paste(colnames(thisLine)[8:9], "_maternal", sep="")
            thisLine$score_paternal <- NA
            thisLine$pvalue_paternal <- NA

            thisLine_gr <- makeGRangesFromDataFrame(thisLine)
            lineIntersect <- as.data.frame(findOverlaps(thisLine_gr,relevantSeqs_gr))
            thisLine$region <- relevantSeqs[lineIntersect[1,2]]

            motif_destruction_info_table <- rbind(motif_destruction_info_table, thisLine)
          }
        }
      }

    }
  }
  return(motif_destruction_info_table)

}


################################################################################
#The following function determines, for each case of motif-disruption, the
#number of reads from the relevant TF mapping to the ancestral and derived alleles.
################################################################################
determine_ancestral_reads <- function(motif_destruction_info_table, motif_destruction_info_table_gr, suppTable1) {
  motif_destruction_info_table_with_ancestral_and_reads <- c()
  weird_entries <- c()
  for (i in 1:nrow(motif_destruction_info_table)) {
    if(i%%100==0) { print(paste(i, nrow(motif_destruction_info_table)))}
    thisEntry <- suppTable1[which(suppTable1$region==motif_destruction_info_table[i,"region"]),]
    thisTF <- as.character(motif_destruction_info_table[i,"motif_alt_id"])
    ################################################################################
    #for each variant, identify the maternal and paternal allele, as well as the
    #ancestral allele and the DAF if the information exists.
    ################################################################################
    theMat_haplo <- unlist(strsplit(thisEntry[1,"hap1"], split=";"))
    theMat_haplo <- matrix(nrow=length(theMat_haplo), ncol=3, data=unlist(strsplit(theMat_haplo, split=",")), byrow=T)

    thePat_haplo <- unlist(strsplit(thisEntry[1,"hap2"], split=";"))
    thePat_haplo <- matrix(nrow=length(thePat_haplo), ncol=3, data=unlist(strsplit(thePat_haplo, split=",")), byrow=T)

    ancestral <- unlist(strsplit(thisEntry[1,"AncAll"], split=";"))
    ancestral <- matrix(nrow=length(ancestral), ncol=2, data=unlist(strsplit(ancestral, split=":", fixed=T)), byrow=T)

    DAF <- unlist(strsplit(thisEntry[1,"DAF"], split=";"))
    DAF <- matrix(nrow=length(DAF), ncol=2, data=unlist(strsplit(DAF, split=":", fixed=T)), byrow=T)

    TF <- strsplit(thisEntry[1,"SigFactors"], split=";")[[1]]
    matReads <- strsplit(thisEntry[1,"hap1_depth"], split=";")[[1]]
    patReads <- strsplit(thisEntry[1,"hap2_depth"], split=";")[[1]]
    pval <- strsplit(thisEntry[1,"pval"], split=";")[[1]]
    TF_mat <- cbind(TF, matReads, patReads, pval)
    TF_mat <- rbind(TF_mat[TF_mat[,1]==thisTF,])

    final_variant_matrix <- cbind(rbind(theMat_haplo[,1:2]), rbind(thePat_haplo[,2:3]))
    final_variant_matrix <- as.data.frame(final_variant_matrix)
    colnames(final_variant_matrix) <- c("position", "mat", "pat", "type")
    final_variant_matrix$Ancestral <- rep(NA, nrow(final_variant_matrix))
    final_variant_matrix$DAF <- rep(NA, nrow(final_variant_matrix))

    for (j in 1:nrow(final_variant_matrix)) {
      if(final_variant_matrix[j,"position"]%in%ancestral[,1]) { final_variant_matrix[j,"Ancestral"] <- ancestral[which(ancestral[,1]==final_variant_matrix[j,"position"]),2]}
      if(final_variant_matrix[j,"position"]%in%DAF[,1]) { final_variant_matrix[j,"DAF"] <- DAF[which(DAF[,1]==final_variant_matrix[j,"position"]),2]}
    }

    final_variant_matrix$TFAffected <- rep("No", (nrow(final_variant_matrix)))
    final_variant_matrix$matReads <- rep(NA, (nrow(final_variant_matrix)))
    final_variant_matrix$patReads <- rep(NA, (nrow(final_variant_matrix)))
    final_variant_matrix$pval <- rep(NA, (nrow(final_variant_matrix)))
    if(nrow(TF_mat)==1) {
      final_variant_matrix$TFAffected <- rep("Yes", nrow(final_variant_matrix))
      final_variant_matrix$matReads <- rep(TF_mat[1,"matReads"], nrow(final_variant_matrix))
      final_variant_matrix$patReads <- rep(TF_mat[1,"patReads"], nrow(final_variant_matrix))
      final_variant_matrix$pval <- rep(TF_mat[1,"pval"], nrow(final_variant_matrix))
    }

    final_variant_matrix$whichAncestral <- rep(NA, nrow(final_variant_matrix))
    final_variant_matrix[which(final_variant_matrix$mat==final_variant_matrix$Ancestral),"whichAncestral"] <- "maternal"
    final_variant_matrix[which(final_variant_matrix$pat==final_variant_matrix$Ancestral),"whichAncestral"] <- "paternal"

    ################################################################################
    #Given this information, identify which variant/variants overlap with the motif
    #of interest create a GR.
    ################################################################################
    final_variant_matrix_gr <- cbind(rep(thisEntry[,1], nrow(final_variant_matrix)), final_variant_matrix$position, final_variant_matrix$position)
    final_variant_matrix_gr <- as.data.frame(final_variant_matrix_gr)
    colnames(final_variant_matrix_gr) <- c("chr", "start", "end")
    final_variant_matrix_gr <- makeGRangesFromDataFrame(final_variant_matrix_gr)

    theIntersect <- as.data.frame(findOverlaps(final_variant_matrix_gr, motif_destruction_info_table_gr[i]))
    if(nrow(theIntersect)>0) {
      for (j in 1:nrow(theIntersect)) {
        thisLine <- cbind(motif_destruction_info_table[i,], final_variant_matrix[theIntersect[j,1],])
        motif_destruction_info_table_with_ancestral_and_reads <- rbind(motif_destruction_info_table_with_ancestral_and_reads, thisLine)
      }
    }
    if(nrow(theIntersect)==0) { weird_entries <- c(weird_entries, i)}

  }
  return(list(motif_destruction_info_table_with_ancestral_and_reads, weird_entries))
}




################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
theMotifs <- args[2]
pval_d1 <- args[3]
pval_d2 <- args[4]
maternal_fasta_d1 <- args[5]
paternal_fasta_d1 <- args[6]
maternal_fasta_d2 <- args[7]
paternal_fasta_d2 <- args[8]
suppTable1 <- args[9]
suppTable2 <- args[10]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#theMotifs <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/jaspar_2023July03/JASPAR2022_CORE_non-redundant_pfms_meme.txt"
#pval_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pval_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#maternal_fasta_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences_maternal.fasta"
#paternal_fasta_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences_paternal.fasta"
#maternal_fasta_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences_maternal.fasta"
#paternal_fasta_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences_paternal.fasta"
#suppTable1 <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/Supplemental_Table_1.txt"
#suppTable2 <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023May26/Supplemental_Table_2.txt"

################################################################################
#Read in motifs
################################################################################

jaspar_motifs <- read_meme(theMotifs)




################################################################################
################################################################################
#Donor 1.
################################################################################
################################################################################

################################################################################
#Read in sequences and significance.
################################################################################

maternal_fasta <- readLines(maternal_fasta_d1)
paternal_fasta <- readLines(paternal_fasta_d1)


pval_d1 <- read.table(pval_d1, header=T, sep="\t", stringsAsFactors=F)
pval_d1 <- pval_d1[-grep("chrX", rownames(pval_d1)),]

sigInput_d1 <- rownames(pval_d1)[which(pval_d1[,"INPUT"]<=0.05)]

pval_d1 <- pval_d1[!rownames(pval_d1)%in%sigInput_d1,]

uniqueTFs <- colnames(pval_d1)
uniqueTFs <- uniqueTFs[order(uniqueTFs)]
uniqueTFs <- uniqueTFs[!uniqueTFs%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "INPUT", "POL2")]

pval_d1 <- pval_d1[,uniqueTFs]

pval_d1 <- pval_d1[which(rowMins(as.matrix(pval_d1))<=0.001),]


info_table_d1 <- c()
for (i in 1:length(uniqueTFs)) {
#for (i in 1:10) {
  alt_jaspar_motifs <- filter_motifs(jaspar_motifs, altname=uniqueTFs[i])
  if(length(alt_jaspar_motifs)>0) {
    print(paste(i, length(uniqueTFs)))

    #Grab sequences of allele-biased binding for this TF.
    #Save them.
    relevantSeqs <- rownames(pval_d1[pval_d1[,uniqueTFs[i]]<=0.001,])
    relevantSeqs <- paste(">", relevantSeqs, sep="")
    rel_grep <- which(maternal_fasta%in%relevantSeqs)
    rel_grep <- c(rel_grep, rel_grep+1)
    rel_grep <- rel_grep[order(rel_grep)]

    this_maternal_fasta <- maternal_fasta[rel_grep]
    this_paternal_fasta <- paternal_fasta[rel_grep]

    saveFile_maternal <- paste(outDir, "Donor1_VarSequences_maternal_TFSig_TFspecific.fasta", sep="")
    write(this_maternal_fasta, saveFile_maternal, sep="\n", append=F)
    saveFile_paternal <- paste(outDir, "Donor1_VarSequences_paternal_TFSig_TFspecific.fasta", sep="")
    write(this_paternal_fasta, saveFile_paternal, sep="\n", append=F)


    #Run fimo on these sequences.
    maternal_Hits <- runFimo(sequences=saveFile_maternal, motifs=alt_jaspar_motifs, meme_path="/cluster/software/meme-5.3.3/bin")
    paternal_Hits <- runFimo(sequences=saveFile_paternal, motifs=alt_jaspar_motifs, meme_path="/cluster/software/meme-5.3.3/bin")

    #Determine motif destruction and disruption.

    relevantSeqs <- gsub("^>", "", relevantSeqs)
    relevantSeqs_matrix <- matrix(nrow=length(relevantSeqs), ncol=3, data=unlist(strsplit(relevantSeqs, split=":|-")), byrow=TRUE)
    relevantSeqs_matrix <- as.data.frame(relevantSeqs_matrix)
    colnames(relevantSeqs_matrix) <- c("chr", "start", "end")
    for (j in 2:3) { relevantSeqs_matrix[,j] <- as.numeric(as.character(relevantSeqs_matrix[,j]))}
    relevantSeqs_gr <- makeGRangesFromDataFrame(relevantSeqs_matrix)
    motif_destruction_info_table <- determine_motif_destruction(maternal_Hits, paternal_Hits, relevantSeqs_gr)
    disruptedNum <- 0
    if(!is.null(nrow(motif_destruction_info_table))) { disruptedNum <- length(unique(motif_destruction_info_table$region))}
    info_table_d1 <- rbind(info_table_d1, c(uniqueTFs[i], length(relevantSeqs), disruptedNum))
  }

}





maternal_fasta <- readLines(maternal_fasta_d2)
paternal_fasta <- readLines(paternal_fasta_d2)


pval_d2 <- read.table(pval_d2, header=T, sep="\t", stringsAsFactors=F)

sigInput_d2 <- rownames(pval_d2)[which(pval_d2[,"INPUT"]<=0.05)]

pval_d2 <- pval_d2[!rownames(pval_d2)%in%sigInput_d2,]

uniqueTFs <- colnames(pval_d2)
uniqueTFs <- uniqueTFs[order(uniqueTFs)]
uniqueTFs <- uniqueTFs[!uniqueTFs%in%c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3", "H3K9AC", "INPUT", "POL2")]

pval_d2 <- pval_d2[,uniqueTFs]

pval_d2 <- pval_d2[which(rowMins(as.matrix(pval_d2))<=0.001),]




info_table_d2 <- c()
for (i in 1:length(uniqueTFs)) {
#for (i in 1:10) {
  alt_jaspar_motifs <- filter_motifs(jaspar_motifs, altname=uniqueTFs[i])
  if(length(alt_jaspar_motifs)>0) {
    print(paste(i, length(uniqueTFs)))

    #Grab sequences of allele-biased binding for this TF.
    #Save them.
    relevantSeqs <- rownames(pval_d2[pval_d2[,uniqueTFs[i]]<=0.001,])
    relevantSeqs <- paste(">", relevantSeqs, sep="")
    rel_grep <- which(maternal_fasta%in%relevantSeqs)
    rel_grep <- c(rel_grep, rel_grep+1)
    rel_grep <- rel_grep[order(rel_grep)]

    this_maternal_fasta <- maternal_fasta[rel_grep]
    this_paternal_fasta <- paternal_fasta[rel_grep]

    saveFile_maternal <- paste(outDir, "Donor1_VarSequences_maternal_TFSig_TFspecific.fasta", sep="")
    write(this_maternal_fasta, saveFile_maternal, sep="\n", append=F)
    saveFile_paternal <- paste(outDir, "Donor1_VarSequences_paternal_TFSig_TFspecific.fasta", sep="")
    write(this_paternal_fasta, saveFile_paternal, sep="\n", append=F)


    #Run fimo on these sequences.
    maternal_Hits <- runFimo(sequences=saveFile_maternal, motifs=alt_jaspar_motifs, meme_path="/cluster/software/meme-5.3.3/bin")
    paternal_Hits <- runFimo(sequences=saveFile_paternal, motifs=alt_jaspar_motifs, meme_path="/cluster/software/meme-5.3.3/bin")

    #Determine motif destruction and disruption.

    relevantSeqs <- gsub("^>", "", relevantSeqs)
    relevantSeqs_matrix <- matrix(nrow=length(relevantSeqs), ncol=3, data=unlist(strsplit(relevantSeqs, split=":|-")), byrow=TRUE)
    relevantSeqs_matrix <- as.data.frame(relevantSeqs_matrix)
    colnames(relevantSeqs_matrix) <- c("chr", "start", "end")
    for (j in 2:3) { relevantSeqs_matrix[,j] <- as.numeric(as.character(relevantSeqs_matrix[,j]))}
    relevantSeqs_gr <- makeGRangesFromDataFrame(relevantSeqs_matrix)
    motif_destruction_info_table <- determine_motif_destruction(maternal_Hits, paternal_Hits, relevantSeqs_gr)
    disruptedNum <- 0
    if(!is.null(nrow(motif_destruction_info_table))) { disruptedNum <- length(unique(motif_destruction_info_table$region))}
    info_table_d2 <- rbind(info_table_d2, c(uniqueTFs[i], length(relevantSeqs), disruptedNum))
  }

}




info_table_d1 <- as.data.frame(info_table_d1)
colnames(info_table_d1) <- c("TF", "ABB_Regions", "Motif_Disrupted_Regions")
for (i in 2:3) {info_table_d1[,i] <- as.numeric(as.character(info_table_d1[,i]))}

info_table_d2 <- as.data.frame(info_table_d2)
colnames(info_table_d2) <- c("TF", "ABB_Regions", "Motif_Disrupted_Regions")
for (i in 2:3) {info_table_d2[,i] <- as.numeric(as.character(info_table_d2[,i]))}


info_table <- merge(info_table_d1, info_table_d2, by="TF")

colnames(info_table) <- c("TF", "Donor1_Biased_Binding", "Donor1_Motif_Disrupted", "Donor2_Biased_Binding", "Donor2_Motif_Disrupted")

info_table$Fraction <- (info_table$Donor1_Motif_Disrupted + info_table$Donor2_Motif_Disrupted) / (info_table$Donor1_Biased_Binding + info_table$Donor2_Biased_Binding)
info_table <- info_table[order(info_table$Fraction, decreasing=T),]

saveFile <- paste(outDir, "Supplemental_Table_4.txt", sep="")
write.table(info_table, saveFile, row.names=F, col.names=T, sep="\t", quote=F)





#
