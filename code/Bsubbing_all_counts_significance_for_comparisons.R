#!/usr/bin/env R
#Bsubbing_all_counts_significance_for_comparisons.R


################################################################################
#This script is used to submit jobs for counting the number of reads supporting
#     one haplotype over another for heterozygous regions and performing binomial
#     tests after mapping has been performed. The script below is tailored to our
#     local compute cluster's specific architecture, and will need to be
#     rewritten for your own compute environment's specific needs if you would
#     like to reproduce its results.
#
#Run under R 4.1.0, but can be used under other versions.
#
#This program takes the following arguments:
# all_hetFiles_dir: Directory in which the sam files with reads mapped to
#     heterozygous regions are stored. It is expected that the file will
#     end with the pattern ".filt.nodup.srt.hetsOnly.sam", and start with
#     a relevant dataset identifier, e.g. "AMY_BCL11A.filt.nodup.srt.hetsOnly.sam"
# allVarSeqs: path to the fasta file containing the heterozygous regions, produced
#     as part of the Building_fasta_variant_sequences.sh script.
#
################################################################################


################################################################################
################################################################################
#Load Libraries and define options.
################################################################################
################################################################################


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
all_hetFiles_dir <- args[1]
allVarSeqs <- args[2]

all_hetFiles <- list.files(all_hetFiles_dir, pattern="filt.nodup.srt.hetsOnly.sam", full.names=T)


for (i in 1:length(all_hetFiles)) {
  hetSam <- all_hetFiles[i]
  theBase <- gsub(".filt.nodup.srt.hetsOnly.sam", "", hetSam)
  logFile <- paste(theBase, "_bsubLog_recount.bsubLog.txt", sep="")
  outputBase <- paste(theBase, "_variantReads_counts", sep="")

  checkFile <- paste(gsub("_counts", "", outputBase), "_binomial_pval.txt", sep="")

  if(!file.exists(checkFile)) {
    theCommand <- paste("bsub -n 3 -R rusage[mem=30000] -We 100:00 -q c7normal -o ", logFile, " -J recounts_vg \"module load cluster/R/4.1.0; Rscript /cluster/home/bmoyers/Scripts/BrainTF_AlleleSpecificBinding/CompilingHaplotypeReads_multipleMappings_vg_v12_2022May16.R ", hetSam, " ", allVarSeqs, " ", outputBase, "\"", sep="")
    system(theCommand)
  }
}
