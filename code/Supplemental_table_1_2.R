#Supplemental_table_1_2.R


################################################################################
#This script is used to produce Supplemental Tables 1 and 2.
#
#Run under R 4.1.0, but can be used under other versions so long as the packages
#     GenomicRanges, matrixStats, and ggplot are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# vcf_saveFile: Path to a VCF-like table that contains all the variants present
#     in either donor, with ancestral allele information added, produced by the
#     script Figure_3B.R .
# allSeqs_d1: path to the fasta file for all variant sequences in donor 1,
#     produced by the script Building_fasta_variant_sequences.sh
# pvals_d1: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# depthTable_mat_d1: A table compiling all TF read  depths for donor 1, haplotype 1,
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depthTable_pat_d1: As depthTable_mat_d1, but for haplotype 2.
# vep_d1_saveFile: a vep-annotated vcf file for donor 1, produced by vep using the
#     vep108.ini file.
# allSeqs_d2: path to the fasta file for all variant sequences in donor 2,
#     produced by the script Building_fasta_variant_sequences.sh
# pvals_d2: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 2, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# depthTable_mat_d2: A table compiling all TF read  depths for donor 2, haplotype 1,
#     when summed across tissues, compiled by the script Summing_depths_across_tissues.R .
# depthTable_pat_d2: As depthTable_mat_d2, but for haplotype 2.
# vep_d2_saveFile: a vep-annotated vcf file for donor 2, produced by vep using the
#     vep108.ini file.
#
################################################################################


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(GenomicRanges)
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
vcf_saveFile <- args[2]
allSeqs_d1 <- args[3]
pvals_d1 <- args[4]
depthTable_mat_d1 <- args[5]
depthTable_pat_d1 <- args[6]
vep_d1_saveFile <- args[7]
allSeqs_d2 <- args[8]
pvals_d2 <- args[9]
depthTable_mat_d2 <- args[10]
depthTable_pat_d2 <- args[11]
vep_d2_saveFile <- args[12]


#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#vcf_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/ensembl/vcf_all_variants_in_donors_1_and_2_with_AncestralAllele_info.txt"
#allSeqs_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#pvals_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#depthTable_mat_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depthTable_pat_d1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#vep_d1_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0001_annotated.tsv.gz"
#allSeqs_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#pvals_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#depthTable_mat_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_maternal_vg_summedAcrossTissues.txt"
#depthTable_pat_d2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypeReadsTable_paternal_vg_summedAcrossTissues.txt"
#vep_d2_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0002_annotated.tsv.gz"





################################################################################
################################################################################
#Donor 1
################################################################################
################################################################################

################################################################################
#Read in a table of all variants.
#Load in p-value table and summed reads tables.
#Identify all cases of min(pval) <=0.05.
#For all cases, identify:
#1) The 2 haplotypes.
#2) All cases for which the variant is significant (excluding INPUT)
#3) The p-value associated with each.
#4) The hap 1 and hap2 reads associated with each.
#5) The overall region considered.
################################################################################

allSeqs <- readLines(allSeqs_d1)
allSeqs <- allSeqs[grep("^>", allSeqs)]
allSeqs <- gsub("^>", "", allSeqs)
allSeqs_Mat <- allSeqs[grep(":M:", allSeqs, fixed=T)]
allSeqs_Pat <- allSeqs[grep(":P:", allSeqs, fixed=T)]
rm(allSeqs)
allSeqs_Mat <- matrix(nrow=length(allSeqs_Mat), ncol=2, data=unlist(strsplit(allSeqs_Mat, split=":M:", fixed=T)), byrow=T)
allSeqs_Pat <- matrix(nrow=length(allSeqs_Pat), ncol=2, data=unlist(strsplit(allSeqs_Pat, split=":P:", fixed=T)), byrow=T)



pvals <- read.table(pvals_d1, header=T, sep="\t", stringsAsFactors=F)

#pvals <- pvals[,colnames(pvals)!="INPUT"]
depthTable_mat <- read.table(depthTable_mat_d1, header=T, sep="\t", stringsAsFactors=F)
depthTable_pat <- read.table(depthTable_pat_d1, header=T, sep="\t", stringsAsFactors=F)

depthTable_mat <- depthTable_mat[,colnames(depthTable_mat)%in%colnames(pvals)]
depthTable_pat <- depthTable_pat[,colnames(depthTable_pat)%in%colnames(pvals)]


toRemove <- grep("chrX", rownames(pvals))
if(length(toRemove)>0) {
  pvals <- pvals[-toRemove,]
  depthTable_mat <- depthTable_mat[-toRemove,]
  depthTable_pat <- depthTable_pat[-toRemove,]
}


theSigs <- which(rowMins(as.matrix(pvals))<=0.05)
pvals <- pvals[theSigs,]
depthTable_mat <- depthTable_mat[theSigs,]
depthTable_pat <- depthTable_pat[theSigs,]

allSeqs_Mat <- allSeqs_Mat[allSeqs_Mat[,1]%in%rownames(pvals),]
allSeqs_Pat <- allSeqs_Pat[allSeqs_Pat[,1]%in%rownames(pvals),]


pvals <- pvals[order(rownames(pvals)),]
depthTable_mat <- depthTable_mat[order(rownames(depthTable_mat)),]
depthTable_pat <- depthTable_pat[order(rownames(depthTable_pat)),]
allSeqs_Mat <- allSeqs_Mat[order(allSeqs_Mat[,1]),]
allSeqs_Pat <- allSeqs_Pat[order(allSeqs_Pat[,1]),]



suppTable_1 <- matrix(nrow=nrow(pvals), ncol=9)
colnames(suppTable_1) <- c("chr", "start", "end", "hap1", "hap2", "SigFactors", "pval", "hap1_depth", "hap2_depth")

for (i in 1:nrow(pvals)) {
  if(i%%1000==0) {print(paste(i, nrow(pvals)))}
  whichSig <- which(as.numeric(pvals[i,])<=0.05)

  theFactors <- colnames(pvals)[whichSig]
  thePvals <- as.numeric(pvals[i,whichSig])
  hap1_depths <- as.numeric(depthTable_mat[i,whichSig])
  hap2_depths <- as.numeric(depthTable_pat[i,whichSig])

  theFactors <- paste(theFactors, collapse=";")
  thePvals <- paste(thePvals, collapse=";")
  hap1_depths <- paste(hap1_depths, collapse=";")
  hap2_depths <- paste(hap2_depths, collapse=";")

  thisHap1 <- allSeqs_Mat[i,2]
  thisHap2 <- allSeqs_Pat[i,2]

  theRegion <- strsplit(rownames(pvals)[i], split=":|-")[[1]]

  thisLine <- c(theRegion, thisHap1, thisHap2, theFactors, thePvals, hap1_depths, hap2_depths)
  thisLine <- rbind(thisLine)
  #thisLine <- as.data.frame(thisLine)
  colnames(thisLine) <- c("chr", "start", "end", "hap1", "hap2", "SigFactors", "pval", "hap1_depth", "hap2_depth")
  #suppTable_1 <- rbind(suppTable_1, thisLine)
  suppTable_1[i,] <- thisLine

}

suppTable_1 <- as.data.frame(suppTable_1)
suppTable_1_gr <- makeGRangesFromDataFrame(suppTable_1)


suppTable_1$hap1 <- gsub(";$", "", suppTable_1$hap1)
suppTable_1$hap2 <- gsub(";$", "", suppTable_1$hap2)




################################################################################
#Determine, for each case, the rsID and ancestral allele, if that information exists.
################################################################################

################################################################################
#First, for every variant, determine if we have ancestral allele  and rsID
#information.
################################################################################

donor_variant_table_d1 <- read.table("/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/phased_variants.vcf.table", header=F, sep="\t", stringsAsFactors=F)
donor_variant_table_d1$name <- paste(donor_variant_table_d1[,1], donor_variant_table_d1[,2], sep="_")

donor_variant_table_d1_gr <- donor_variant_table_d1[,c(1,2,2)]
colnames(donor_variant_table_d1_gr) <- c("chr", "start", "end")
donor_variant_table_d1_gr <- as.data.frame(donor_variant_table_d1_gr)
donor_variant_table_d1_gr <- makeGRangesFromDataFrame(donor_variant_table_d1_gr)

theIntersect <- as.data.frame(findOverlaps(donor_variant_table_d1_gr, suppTable_1_gr))

donor_variant_table_d1 <- donor_variant_table_d1[unique(theIntersect[,1]),]

donor_variant_table_d1_gr <- donor_variant_table_d1[,c(1,2,2)]
colnames(donor_variant_table_d1_gr) <- c("chr", "start", "end")
donor_variant_table_d1_gr <- as.data.frame(donor_variant_table_d1_gr)
donor_variant_table_d1_gr <- makeGRangesFromDataFrame(donor_variant_table_d1_gr)


vcf_in_our_donors_with_ancestral_allele <- read.table(vcf_saveFile, header=T, sep="\t", stringsAsFactors=F)


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

vcf_in_our_donors_with_ancestral_allele_df <- vcf_in_our_donors_with_ancestral_allele[,c(1,2,2)]
colnames(vcf_in_our_donors_with_ancestral_allele_df) <- c("chr", "start", "end")
vcf_in_our_donors_with_ancestral_allele_gr <- makeGRangesFromDataFrame(vcf_in_our_donors_with_ancestral_allele_df)



mini_vcf <- vcf_in_our_donors_with_ancestral_allele[vcf_in_our_donors_with_ancestral_allele$name%in%donor_variant_table_d1$name,]


mini_vcf <- mini_vcf[order(mini_vcf$name),]
donor_variant_table_d1 <- donor_variant_table_d1[order(donor_variant_table_d1$name),]


these_rsID <- rep(NA, nrow(donor_variant_table_d1))
these_ancAll <- rep(NA, nrow(donor_variant_table_d1))

currentRow <- 1
for (i in 1:nrow(mini_vcf)) {
  if(i%%100==0) {print(paste(i, nrow(mini_vcf)))}
  while(mini_vcf[i,"name"]!=donor_variant_table_d1[currentRow,"name"]) {
    currentRow <- currentRow+1
    if(currentRow>nrow(donor_variant_table_d1)) {
      print("Ran Out of Rows!")
      break()
    }
  }

  donor_variant_alleles <- unlist(strsplit(as.character(donor_variant_table_d1[currentRow,4:5]), split=","))
  vcf_variant_alleles <- unlist(strsplit(as.character(mini_vcf[i,"ALT"]), split=","))
  vcf_reference_alleles <- unlist(strsplit(as.character(mini_vcf[i,"REF"]), split=","))

  if(length(donor_variant_alleles[donor_variant_alleles%in%vcf_variant_alleles])>0 && length(donor_variant_alleles[donor_variant_alleles%in%c(vcf_variant_alleles, vcf_reference_alleles)])==length(donor_variant_alleles)) {
    if(!is.na(these_rsID[currentRow])) { these_rsID[currentRow] <- paste(these_rsID[currentRow], mini_vcf[i,"ID"], sep=",")}
    if(is.na(these_rsID[currentRow])) { these_rsID[currentRow] <- mini_vcf[i,"ID"] }

    if(!is.na(these_ancAll[currentRow])) {
      if(these_ancAll[currentRow]!=mini_vcf[i,"AncAll"]) {
        these_ancAll[currentRow] <- paste(these_ancAll[currentRow], mini_vcf[i,"AncAll"], sep=",")
      }
    }

    if(is.na(these_ancAll[currentRow])) { these_ancAll[currentRow] <- mini_vcf[i,"AncAll"] }
  }

}


donor_variant_table_d1$rsID <- these_rsID
donor_variant_table_d1$AncAll <- these_ancAll

donor_variant_table_d1_gr <- donor_variant_table_d1[,c(1,2,2)]
colnames(donor_variant_table_d1_gr) <- c("chr", "start", "end")
donor_variant_table_d1_gr <- as.data.frame(donor_variant_table_d1_gr)
donor_variant_table_d1_gr <- makeGRangesFromDataFrame(donor_variant_table_d1_gr)


################################################################################
#Next, we need to determine the derived allele frequency at this location.
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


#


name1 <- rep(NA, nrow(donor_variant_table_d1))
name2 <- rep(NA, nrow(donor_variant_table_d1))
for (i in 1:nrow(donor_variant_table_d1)) {
  if(i%%1000==0) {print(paste(i, nrow(donor_variant_table_d1)))}
  twoOptions <- donor_variant_table_d1[i,c(4, 5)]
  theAncestralAllele <- donor_variant_table_d1[i,"AncAll"]
  theDerivedAllele <- twoOptions[twoOptions!=theAncestralAllele]
  if(length(theDerivedAllele)==1) {
    this_name1 <- paste(donor_variant_table_d1[i,1], "_", donor_variant_table_d1[i,2], "_", theAncestralAllele, "/", theDerivedAllele, sep="")
    name1[i] <- this_name1
    this_name2 <- paste(donor_variant_table_d1[i,1], "_", donor_variant_table_d1[i,2], "_", theDerivedAllele, "/", theAncestralAllele, sep="")
    name2[i] <- this_name2
  }

}


donor_variant_table_d1$name1 <- name1
donor_variant_table_d1$name2 <- name2


match1 <- match(donor_variant_table_d1$name1, vep_d1$Uploaded_variation)
match2 <- match(donor_variant_table_d1$name2, vep_d1$Uploaded_variation)


DAF1 <- rep(NA, nrow(donor_variant_table_d1))
DAF1 <- vep_d1[match1, "gnomad3_AF"]

DAF2 <- rep(NA, nrow(donor_variant_table_d1))
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


donor_variant_table_d1$DAF <- finalDAF




#donor_variant_table_d1$name1 <- paste(donor_variant_table_d1$V1, "_", donor_variant_table_d1$V2, "_", donor_variant_table_d1$V4, "/", donor_variant_table_d1$V5, sep="")
#match1 <- match(donor_variant_table_d1$name1, vep_d1$Uploaded_variation)
#DAF <- rep(NA, nrow(donor_variant_table_d1))
#DAF <- vep_d1[match1, "gnomad3_AF"]
#donor_variant_table_d1$DAF <- DAF



################################################################################
#Given that information, intersect the donor variant table with the suppTable1.
#For each intersect, record for the specific location ancestral alleles and rsIDs.
################################################################################


theIntersect <- as.data.frame(findOverlaps(donor_variant_table_d1_gr, suppTable_1_gr))


final_rsIDs <- rep(NA, nrow(suppTable_1))
final_ancAll <- rep(NA, nrow(suppTable_1))
final_DAF <- rep(NA, nrow(suppTable_1))

for (i in 1:nrow(theIntersect)) {
  if(i%%100==0) {print(paste(i, nrow(theIntersect)))}
  thisSupp <- theIntersect[i,2]
  thisVCF <- theIntersect[i,1]

  if(!is.na(donor_variant_table_d1[thisVCF,"rsID"])) {
    if(!is.na(final_rsIDs[thisSupp])) { final_rsIDs[thisSupp] <- paste(final_rsIDs[thisSupp], paste(donor_variant_table_d1[thisVCF,2], donor_variant_table_d1[thisVCF,"rsID"], sep=":"), sep=";")}
    if(is.na(final_rsIDs[thisSupp])) { final_rsIDs[thisSupp] <- paste(donor_variant_table_d1[thisVCF,2], donor_variant_table_d1[thisVCF,"rsID"], sep=":")}
  }

  if(!is.na(donor_variant_table_d1[thisVCF,"DAF"])) {
    if(!is.na(final_DAF[thisSupp])) { final_DAF[thisSupp] <- paste(final_DAF[thisSupp], paste(donor_variant_table_d1[thisVCF,2], donor_variant_table_d1[thisVCF,"DAF"], sep=":"), sep=";")}
    if(is.na(final_DAF[thisSupp])) { final_DAF[thisSupp] <- paste(donor_variant_table_d1[thisVCF,2], donor_variant_table_d1[thisVCF,"DAF"], sep=":")}
  }

  if(!is.na(donor_variant_table_d1[thisVCF,"AncAll"])) {
    if(!is.na(final_ancAll[thisSupp])) {
      theAdd <- paste(donor_variant_table_d1[thisVCF,2], donor_variant_table_d1[thisVCF,"AncAll"], sep=":")
      theCompare <- unlist(strsplit(final_ancAll[thisSupp], split=";"))
      if(!theAdd%in%theCompare) { final_ancAll[thisSupp] <- paste(final_ancAll[thisSupp], theAdd, sep=";") }
    }
    if(is.na(final_ancAll[thisSupp])) { final_ancAll[thisSupp] <- paste(donor_variant_table_d1[thisVCF,2], donor_variant_table_d1[thisVCF,"AncAll"], sep=":")}
  }

}


suppTable_1$rsID <- final_rsIDs
suppTable_1$AncAll <- final_ancAll
suppTable_1$DAF <- final_DAF




suppTable_1 <- suppTable_1[,c("chr", "start", "end", "hap1", "hap2", "rsID", "AncAll", "DAF", "SigFactors", "pval", "hap1_depth", "hap2_depth")]



################################################################################
#Save the supplementary table
################################################################################


saveFile <- paste(outDir, "Supplemental_Table_1.txt", sep="")
write.table(suppTable_1, saveFile, row.names=F, col.names=T, sep="\t", quote=F)







################################################################################
################################################################################
#Donor 2
################################################################################
################################################################################

################################################################################
#Read in a table of all variants.
#Load in p-value table and summed reads tables.
#Identify all cases of min(pval) <=0.05.
#For all cases, identify:
#1) The 2 haplotypes.
#2) All cases for which the variant is significant (excluding INPUT)
#3) The p-value associated with each.
#4) The hap 1 and hap2 reads associated with each.
#5) The overall region considered.
################################################################################

allSeqs <- readLines(allSeqs_d2)
allSeqs <- allSeqs[grep("^>", allSeqs)]
allSeqs <- gsub("^>", "", allSeqs)
allSeqs_Mat <- allSeqs[grep(":M:", allSeqs, fixed=T)]
allSeqs_Pat <- allSeqs[grep(":P:", allSeqs, fixed=T)]
rm(allSeqs)
allSeqs_Mat <- matrix(nrow=length(allSeqs_Mat), ncol=2, data=unlist(strsplit(allSeqs_Mat, split=":M:", fixed=T)), byrow=T)
allSeqs_Pat <- matrix(nrow=length(allSeqs_Pat), ncol=2, data=unlist(strsplit(allSeqs_Pat, split=":P:", fixed=T)), byrow=T)



pvals <- read.table(pvals_d2, header=T, sep="\t", stringsAsFactors=F)

#pvals <- pvals[,colnames(pvals)!="INPUT"]

depthTable_mat <- read.table(depthTable_mat_d2, header=T, sep="\t", stringsAsFactors=F)
depthTable_pat <- read.table(depthTable_pat_d2, header=T, sep="\t", stringsAsFactors=F)

depthTable_mat <- depthTable_mat[,colnames(depthTable_mat)%in%colnames(pvals)]
depthTable_pat <- depthTable_pat[,colnames(depthTable_pat)%in%colnames(pvals)]


toRemove <- grep("chrX", rownames(pvals))
if(length(toRemove)>0) {
  pvals <- pvals[-toRemove,]
  depthTable_mat <- depthTable_mat[-toRemove,]
  depthTable_pat <- depthTable_pat[-toRemove,]
}


theSigs <- which(rowMins(as.matrix(pvals))<=0.05)
pvals <- pvals[theSigs,]
depthTable_mat <- depthTable_mat[theSigs,]
depthTable_pat <- depthTable_pat[theSigs,]

allSeqs_Mat <- allSeqs_Mat[allSeqs_Mat[,1]%in%rownames(pvals),]
allSeqs_Pat <- allSeqs_Pat[allSeqs_Pat[,1]%in%rownames(pvals),]


pvals <- pvals[order(rownames(pvals)),]
depthTable_mat <- depthTable_mat[order(rownames(depthTable_mat)),]
depthTable_pat <- depthTable_pat[order(rownames(depthTable_pat)),]
allSeqs_Mat <- allSeqs_Mat[order(allSeqs_Mat[,1]),]
allSeqs_Pat <- allSeqs_Pat[order(allSeqs_Pat[,1]),]



suppTable_2 <- matrix(nrow=nrow(pvals), ncol=9)
colnames(suppTable_2) <- c("chr", "start", "end", "hap1", "hap2", "SigFactors", "pval", "hap1_depth", "hap2_depth")

for (i in 1:nrow(pvals)) {
  if(i%%1000==0) {print(paste(i, nrow(pvals)))}
  whichSig <- which(as.numeric(pvals[i,])<=0.05)

  theFactors <- colnames(pvals)[whichSig]
  thePvals <- as.numeric(pvals[i,whichSig])
  hap1_depths <- as.numeric(depthTable_mat[i,whichSig])
  hap2_depths <- as.numeric(depthTable_pat[i,whichSig])

  theFactors <- paste(theFactors, collapse=";")
  thePvals <- paste(thePvals, collapse=";")
  hap1_depths <- paste(hap1_depths, collapse=";")
  hap2_depths <- paste(hap2_depths, collapse=";")

  thisHap1 <- allSeqs_Mat[i,2]
  thisHap2 <- allSeqs_Pat[i,2]

  theRegion <- strsplit(rownames(pvals)[i], split=":|-")[[1]]

  thisLine <- c(theRegion, thisHap1, thisHap2, theFactors, thePvals, hap1_depths, hap2_depths)
  thisLine <- rbind(thisLine)
  #thisLine <- as.data.frame(thisLine)
  colnames(thisLine) <- c("chr", "start", "end", "hap1", "hap2", "SigFactors", "pval", "hap1_depth", "hap2_depth")
  suppTable_2[i,] <- thisLine

}

suppTable_2 <- as.data.frame(suppTable_2)
suppTable_2_gr <- makeGRangesFromDataFrame(suppTable_2)


suppTable_2$hap1 <- gsub(";$", "", suppTable_2$hap1)
suppTable_2$hap2 <- gsub(";$", "", suppTable_2$hap2)




################################################################################
#Determine, for each case, if a snpID exists (VEP file? VCF file?)
#Get the ancestral allele.  Get Derived allele frequency.  So let's do VEP.
################################################################################

################################################################################
#First, for every variant, determine if we have ancestral allele  and rsID
#information.
################################################################################

donor_variant_table_d2 <- read.table("/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/phased_variants.vcf.table", header=F, sep="\t", stringsAsFactors=F)
donor_variant_table_d2$name <- paste(donor_variant_table_d2[,1], donor_variant_table_d2[,2], sep="_")

donor_variant_table_d2_gr <- donor_variant_table_d2[,c(1,2,2)]
colnames(donor_variant_table_d2_gr) <- c("chr", "start", "end")
donor_variant_table_d2_gr <- as.data.frame(donor_variant_table_d2_gr)
donor_variant_table_d2_gr <- makeGRangesFromDataFrame(donor_variant_table_d2_gr)

theIntersect <- as.data.frame(findOverlaps(donor_variant_table_d2_gr, suppTable_2_gr))

donor_variant_table_d2 <- donor_variant_table_d2[unique(theIntersect[,1]),]

donor_variant_table_d2_gr <- donor_variant_table_d2[,c(1,2,2)]
colnames(donor_variant_table_d2_gr) <- c("chr", "start", "end")
donor_variant_table_d2_gr <- as.data.frame(donor_variant_table_d2_gr)
donor_variant_table_d2_gr <- makeGRangesFromDataFrame(donor_variant_table_d2_gr)


vcf_in_our_donors_with_ancestral_allele <- read.table(vcf_saveFile, header=T, sep="\t", stringsAsFactors=F)


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

vcf_in_our_donors_with_ancestral_allele_df <- vcf_in_our_donors_with_ancestral_allele[,c(1,2,2)]
colnames(vcf_in_our_donors_with_ancestral_allele_df) <- c("chr", "start", "end")
vcf_in_our_donors_with_ancestral_allele_gr <- makeGRangesFromDataFrame(vcf_in_our_donors_with_ancestral_allele_df)



mini_vcf <- vcf_in_our_donors_with_ancestral_allele[vcf_in_our_donors_with_ancestral_allele$name%in%donor_variant_table_d2$name,]


mini_vcf <- mini_vcf[order(mini_vcf$name),]
donor_variant_table_d2 <- donor_variant_table_d2[order(donor_variant_table_d2$name),]


these_rsID <- rep(NA, nrow(donor_variant_table_d2))
these_ancAll <- rep(NA, nrow(donor_variant_table_d2))

currentRow <- 1
for (i in 1:nrow(mini_vcf)) {
  if(i%%100==0) {print(paste(i, nrow(mini_vcf)))}
  while(mini_vcf[i,"name"]!=donor_variant_table_d2[currentRow,"name"]) {
    currentRow <- currentRow+1
    if(currentRow>nrow(donor_variant_table_d2)) {
      print("Ran Out of Rows!")
      break()
    }
  }

  donor_variant_alleles <- unlist(strsplit(as.character(donor_variant_table_d2[currentRow,4:5]), split=","))
  vcf_variant_alleles <- unlist(strsplit(as.character(mini_vcf[i,"ALT"]), split=","))
  vcf_reference_alleles <- unlist(strsplit(as.character(mini_vcf[i,"REF"]), split=","))

  if(length(donor_variant_alleles[donor_variant_alleles%in%vcf_variant_alleles])>0 && length(donor_variant_alleles[donor_variant_alleles%in%c(vcf_variant_alleles, vcf_reference_alleles)])==length(donor_variant_alleles)) {
    if(!is.na(these_rsID[currentRow])) { these_rsID[currentRow] <- paste(these_rsID[currentRow], mini_vcf[i,"ID"], sep=",")}
    if(is.na(these_rsID[currentRow])) { these_rsID[currentRow] <- mini_vcf[i,"ID"] }

    if(!is.na(these_ancAll[currentRow])) {
      if(these_ancAll[currentRow]!=mini_vcf[i,"AncAll"]) {
        these_ancAll[currentRow] <- paste(these_ancAll[currentRow], mini_vcf[i,"AncAll"], sep=",")
      }
    }

    if(is.na(these_ancAll[currentRow])) { these_ancAll[currentRow] <- mini_vcf[i,"AncAll"] }
  }


}


donor_variant_table_d2$rsID <- these_rsID
donor_variant_table_d2$AncAll <- these_ancAll

donor_variant_table_d2_gr <- donor_variant_table_d2[,c(1,2,2)]
colnames(donor_variant_table_d2_gr) <- c("chr", "start", "end")
donor_variant_table_d2_gr <- as.data.frame(donor_variant_table_d2_gr)
donor_variant_table_d2_gr <- makeGRangesFromDataFrame(donor_variant_table_d2_gr)


################################################################################
#Next, we need to determine the derived allele frequency at this location.
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


#


name1 <- rep(NA, nrow(donor_variant_table_d2))
name2 <- rep(NA, nrow(donor_variant_table_d2))
for (i in 1:nrow(donor_variant_table_d2)) {
  if(i%%1000==0) {print(paste(i, nrow(donor_variant_table_d2)))}
  twoOptions <- donor_variant_table_d2[i,c(4, 5)]
  theAncestralAllele <- donor_variant_table_d2[i,"AncAll"]
  theDerivedAllele <- twoOptions[twoOptions!=theAncestralAllele]
  if(length(theDerivedAllele)==1) {
    this_name1 <- paste(donor_variant_table_d2[i,1], "_", donor_variant_table_d2[i,2], "_", theAncestralAllele, "/", theDerivedAllele, sep="")
    name1[i] <- this_name1
    this_name2 <- paste(donor_variant_table_d2[i,1], "_", donor_variant_table_d2[i,2], "_", theDerivedAllele, "/", theAncestralAllele, sep="")
    name2[i] <- this_name2
  }

}


donor_variant_table_d2$name1 <- name1
donor_variant_table_d2$name2 <- name2


match1 <- match(donor_variant_table_d2$name1, vep_d2$Uploaded_variation)
match2 <- match(donor_variant_table_d2$name2, vep_d2$Uploaded_variation)


DAF1 <- rep(NA, nrow(donor_variant_table_d2))
DAF1 <- vep_d2[match1, "gnomad3_AF"]

DAF2 <- rep(NA, nrow(donor_variant_table_d2))
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


donor_variant_table_d2$DAF <- finalDAF









################################################################################
#Given that information, intersect the donor variant table with the suppTable1.
#For each intersect, record for the specific location ancestral alleles and rsIDs.
################################################################################


theIntersect <- as.data.frame(findOverlaps(donor_variant_table_d2_gr, suppTable_2_gr))


final_rsIDs <- rep(NA, nrow(suppTable_2))
final_ancAll <- rep(NA, nrow(suppTable_2))
final_DAF <- rep(NA, nrow(suppTable_2))

for (i in 1:nrow(theIntersect)) {
  if(i%%100==0) {print(paste(i, nrow(theIntersect)))}
  thisSupp <- theIntersect[i,2]
  thisVCF <- theIntersect[i,1]

  if(!is.na(donor_variant_table_d2[thisVCF,"rsID"])) {
    if(!is.na(final_rsIDs[thisSupp])) { final_rsIDs[thisSupp] <- paste(final_rsIDs[thisSupp], paste(donor_variant_table_d2[thisVCF,2], donor_variant_table_d2[thisVCF,"rsID"], sep=":"), sep=";")}
    if(is.na(final_rsIDs[thisSupp])) { final_rsIDs[thisSupp] <- paste(donor_variant_table_d2[thisVCF,2], donor_variant_table_d2[thisVCF,"rsID"], sep=":")}
  }

  if(!is.na(donor_variant_table_d2[thisVCF,"DAF"])) {
    if(!is.na(final_DAF[thisSupp])) { final_DAF[thisSupp] <- paste(final_DAF[thisSupp], paste(donor_variant_table_d2[thisVCF,2], donor_variant_table_d2[thisVCF,"DAF"], sep=":"), sep=";")}
    if(is.na(final_DAF[thisSupp])) { final_DAF[thisSupp] <- paste(donor_variant_table_d2[thisVCF,2], donor_variant_table_d2[thisVCF,"DAF"], sep=":")}
  }

  if(!is.na(donor_variant_table_d2[thisVCF,"AncAll"])) {
    if(!is.na(final_ancAll[thisSupp])) {
      theAdd <- paste(donor_variant_table_d2[thisVCF,2], donor_variant_table_d2[thisVCF,"AncAll"], sep=":")
      theCompare <- unlist(strsplit(final_ancAll[thisSupp], split=";"))
      if(!theAdd%in%theCompare) { final_ancAll[thisSupp] <- paste(final_ancAll[thisSupp], theAdd, sep=";") }
    }
    if(is.na(final_ancAll[thisSupp])) { final_ancAll[thisSupp] <- paste(donor_variant_table_d2[thisVCF,2], donor_variant_table_d2[thisVCF,"AncAll"], sep=":")}
  }

}


suppTable_2$rsID <- final_rsIDs
suppTable_2$AncAll <- final_ancAll
suppTable_2$DAF <- final_DAF




suppTable_2 <- suppTable_2[,c("chr", "start", "end", "hap1", "hap2", "rsID", "AncAll", "DAF", "SigFactors", "pval", "hap1_depth", "hap2_depth")]



################################################################################
#Save the supplementary table
################################################################################


saveFile <- paste(outDir, "Supplemental_Table_2.txt", sep="")
write.table(suppTable_2, saveFile, row.names=F, col.names=T, sep="\t", quote=F)







#








#
