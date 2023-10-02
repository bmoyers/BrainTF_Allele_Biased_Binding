#!/usr/bin/env R
#VCF_to_VCFTable.R


################################################################################
#This script is used to convert a phased VCF to a table that can more easily be
#parsed and used by a later script, Construct_fasta_variant_sequences.R
#
#Run under R 3.4.3, but can be used under other versions.  No libraries required.
#
#This program takes the following arguments:
# VCF: Path to the phased VCF to be parsed, Provided via the following DOI:
#     https://doi.org/10.7303/syn4921369
# VCF_TABLE: Path to which the parsed VCF table will be written.
#
################################################################################


################################################################################
################################################################################
#Load Libraries and define options.
################################################################################
################################################################################

options("scipen"=100)


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
VCF <- args[1]
VCF_TABLE <- args[2]

con <- file(VCF, open="r")

#Reading in the whole VCF takes a while, and is resource-intensive.  Instead,
#I read line-by-line.  As long as there's still another line, read it in and check
#that the variant passed quality filters, pick out relevant information about the
#variant and store it in a table.

current.line <- 1
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(length(grep("^#", line))==0) {
    line <- strsplit(line, split="\t")[[1]]
    if(line[7]=="PASS") {
      thisFormat <- strsplit(line[9], split=":", fixed=T)[[1]]
      thisInfo <- strsplit(line[10], split=":", fixed=T)[[1]]
      thisPhaseSet <- thisInfo[which(thisFormat=="PS")]
      thisGenotype <- thisInfo[which(thisFormat=="GT")]
      finalLine <- paste(c(line[1:5], thisPhaseSet, thisGenotype), collapse="\t")
      write(finalLine, VCF_TABLE, append=T, sep="\n")
    }
  }
  current.line <- current.line+1
}
