#Figure_1B_1C_Supplemental_1_2_3.R


################################################################################
#This script is used to produce Figures 1B, 1C, and Supplemental Figures 1, 2, and 3.
#
#Run under R 4.1.0, but can be used under other versions so long as the packages
#     GenomicRanges, and ggplot2, ComplexHeatmap, and circlize are available.
#
#This program takes the following arguments:
# outDir: Directory to which figures should be saved.
# vg_directory: Directory to which the vg mappings and counts quantifications
#     were written.
# bowtie2_directory: Directory to which the bowtie mappings, pileups, and count
#     quantifications were written.
# variants_rep1: Fasta file containing all variant sequences for donor 1,
#     produced by the script Building_fasta_variant_sequences.sh
# variants_rep2: Fasta file containing all variant sequences for donor 2,
#     produced by the script Building_fasta_variant_sequences.sh
# pval_table_1: A table compiling all p-values for each TF over all variants
#     when reads are summed across tissues for donor 1, produced by the script
#     CompilingHaplotypePvalsTable_CheckNames_summedAcrossTissues.R
# pval_table_2: As pval_table_1, but for donor 2.
#
#
################################################################################


#module load cluster/R/4.1.0
#cd /cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/ComparingBryanSteph


################################################################################
################################################################################
#Load Libraries
################################################################################
################################################################################

library(ggplot2)
library(GenomicRanges)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)

################################################################################
################################################################################
#Define Functions
################################################################################
################################################################################


################################################################################
#This function pulls out the relevant experiment names.
################################################################################
find_all_experiment_names <- function(bowtie2_results) {
  theExprs <- c()
  for (i in 1:length(bowtie2_results)) {
    thisResults <- strsplit(bowtie2_results[i], split="/")[[1]]
    thisResults <- thisResults[length(thisResults)]
    thisResults <- strsplit(thisResults, split="_variantReads_")[[1]][1]
    theExprs <- c(theExprs, thisResults)
  }
  if(length(grep("Control", theExprs))>0) { theExprs <- theExprs[-grep("Control", theExprs)]}
  return(theExprs)
}

################################################################################
#This function determines the statistical significance of each variant in the
#two methods.
################################################################################
getSignificances <- function(this_vg_set, this_bowtie2) {
  belleExpr_mat <- read.table(this_vg_set[1], header=T, sep="\t", stringsAsFactors=F)
  belleExpr_pat <- read.table(this_vg_set[2], header=T, sep="\t", stringsAsFactors=F)

  bowtieExpr_mat <- read.table(this_bowtie2[1], header=T, sep="\t", stringsAsFactors=F)
  bowtieExpr_pat <- read.table(this_bowtie2[2], header=T, sep="\t", stringsAsFactors=F)


  thePvals <- rep(1, nrow(belleExpr_mat))
  theDepths <- cbind(rowSums(cbind(belleExpr_mat, belleExpr_pat)))
  test_vector <- floor(apply(cbind(belleExpr_mat, belleExpr_pat), 1, max))
  to_replace <- which(theDepths[,1]>0)
  bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
  thePvals_replace <- mapply(bt, test_vector[to_replace], theDepths[to_replace,1])
  thePvals[to_replace] <- thePvals_replace

  thePvals_default <- rep(1, nrow(bowtieExpr_mat))
  theDepths_default <- cbind(rowSums(cbind(bowtieExpr_mat, bowtieExpr_pat)))
  test_vector <- floor(apply(cbind(bowtieExpr_mat, bowtieExpr_pat), 1, max))
  to_replace <- which(theDepths_default[,1]>0)
  bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
  thePvals_default_replace <- mapply(bt, test_vector[to_replace], theDepths_default[to_replace,1])
  thePvals_default[to_replace] <- thePvals_default_replace

  finalMat <- matrix(nrow=length(thePvals), ncol=3, data=unlist(strsplit(rownames(belleExpr_mat), split=":|-")), byrow=T)
  finalMat <- cbind(finalMat, belleExpr_mat[,1], belleExpr_pat[,1], theDepths, thePvals, bowtieExpr_mat[,1], bowtieExpr_pat[,1], theDepths_default, thePvals_default)
  finalMat <- as.data.frame(finalMat)
  finalMat[,1] <- as.character(finalMat[,1])
  finalMat[,2] <- as.numeric(as.character(finalMat[,2]))
  finalMat[,3] <- as.numeric(as.character(finalMat[,3]))
  finalMat[,4] <- as.numeric(as.character(finalMat[,4]))
  finalMat[,5] <- as.numeric(as.character(finalMat[,5]))
  finalMat[,6] <- as.numeric(as.character(finalMat[,6]))
  finalMat[,7] <- as.numeric(as.character(finalMat[,7]))
  finalMat[,8] <- as.numeric(as.character(finalMat[,8]))
  finalMat[,9] <- as.numeric(as.character(finalMat[,9]))
  finalMat[,10] <- as.numeric(as.character(finalMat[,10]))
  finalMat[,11] <- as.numeric(as.character(finalMat[,11]))
  colnames(finalMat) <- c("chr", "start", "end", "matCount", "patCount", "totalDepth", "pval", "def_matCount", "def_patCount", "def_totalDepth", "default_pval")

  #firstSig <- which(finalMat[,"pval"]<0.05)
  #secondSig <- which(finalMat[,"default_pval"]<0.05)
  #allSig <- unique(c(firstSig, secondSig))
  #allSig <- allSig[order(allSig)]
  #finalMat <- finalMat[allSig,]
  firstDepth <- which(theDepths>=6)
  secondDepth <- which(theDepths_default>=6)
  allDepth <- unique(c(firstDepth, secondDepth))
  allDepth <- allDepth[order(allDepth)]
  finalMat <- finalMat[allDepth,]
  return(finalMat)
}


################################################################################
#This function grabs a quick summary of a variant for each of the two
#experiments, whether they were significant in the experiments or not, etc.
################################################################################
get_variant_descriptions <- function(this_vg_bowtie, thisExpr, variants_rep1, variants_rep2) {
  this_experiment <- strsplit(thisExpr, split="_")[[1]]
  this_rep <- this_experiment[3]
  if(this_rep=="rep1") { thisVars <- variants_rep1 }
  if(this_rep=="rep2") { thisVars <- variants_rep2 }

  region_names <- paste(this_vg_bowtie[,1], ":", this_vg_bowtie[,2], "-", this_vg_bowtie[,3], sep="")
  var_names <- paste(thisVars[,1], ":", thisVars[,2], "-", thisVars[,3], sep="")

  theMatches <- match(region_names, var_names)
  var_descrptions <- thisVars[theMatches,4:5]

  finalMat <- cbind(this_vg_bowtie, var_descrptions)

}

################################################################################
#this function parses the variant structure identifying the location of each
#variant as well as identifying which was the maternal and paternal allele in
#each case.
################################################################################
parseVariants <- function(theFasta) {
  theVars <- readLines(theFasta)
  theVars <- gsub(">", "", theVars)
  theMats <- theVars[grep(":M:", theVars, fixed=T)]
  thePats <- theVars[grep(":P:", theVars, fixed=T)]

  theMats <- matrix(nrow=length(theMats), ncol=2, data=unlist(strsplit(theMats, split=":M:", fixed=T)), byrow=T)
  thePats <- matrix(nrow=length(thePats), ncol=2, data=unlist(strsplit(thePats, split=":P:", fixed=T)), byrow=T)
  theLocs <- matrix(nrow=nrow(theMats), ncol=3, data=unlist(strsplit(theMats[,1], split=":|-")), byrow=T)

  finalMat <- as.data.frame(cbind(theLocs, theMats[,2], thePats[,2]))
  finalMat[,1] <- as.character(finalMat[,1])
  finalMat[,2] <- as.numeric(as.character(finalMat[,2]))
  finalMat[,3] <- as.numeric(as.character(finalMat[,3]))
  finalMat[,4] <- as.character(finalMat[,4])
  finalMat[,5] <- as.character(finalMat[,5])
  colnames(finalMat) <- c("chr", "start", "end", "matVars", "patVars")

  return(finalMat)
}

################################################################################
#this function identifies, in each case, which variant is the maternal and
#which is the paternal.  It then determines the reference allele frequency
#for both vg and bowtie2.  Note that this can only be done for variants
#which either A) have a single location of change or B) the reference allele
#occurs on the same haplotype for all variants in the block.
################################################################################
get_refAlt_info <- function(this_vg_bowtie, thisExpr, thePileups) {
  #thisExpr <- all_experiments[1]
  #this_vg_bowtie_safe <- this_vg_bowtie
  thisPileup <- thePileups[grep(thisExpr, thePileups)]
  thisPileup <- read.table(thisPileup, header=F, sep="\t", stringsAsFactors=F)
  colnames(thisPileup) <- c("chr", "pos", "NA1", "Ref", "NA2", "Depth", "NA3", "Info", "GT", "Hap")
  thisPileup$name <- paste(thisPileup$chr, thisPileup$pos, thisPileup$Ref, sep="_")

  ref_haplotype <- rep(NA, nrow(this_vg_bowtie))
  ref_haplotype_strict <- rep(NA, nrow(this_vg_bowtie))
  for (k in 1:nrow(this_vg_bowtie)) {
    #print(paste(k, nrow(this_vg_bowtie)))
    if(k%%1000==0) { print(paste(k, nrow(this_vg_bowtie)))}
    thisSet_mat <- strsplit(this_vg_bowtie[k,"matVars"], split=";")[[1]]
    thisSet_pat <- strsplit(this_vg_bowtie[k,"patVars"], split=";")[[1]]

    thisSet_mat <- thisSet_mat[grep("HET", thisSet_mat)]
    thisSet_pat <- thisSet_pat[grep("HET", thisSet_pat)]


    for (l in 1:length(thisSet_mat)) {
      thisSet_mat[l] <- paste(this_vg_bowtie[k,"chr"], paste(strsplit(thisSet_mat[l], split=",")[[1]][1:2], collapse="_"), sep="_")
      thisSet_pat[l] <- paste(this_vg_bowtie[k,"chr"], paste(strsplit(thisSet_pat[l], split=",")[[1]][1:2], collapse="_"), sep="_")
    }

    a <- length(thisSet_mat[thisSet_mat%in%thisPileup$name])
    b <- length(thisSet_pat[thisSet_pat%in%thisPileup$name])
    if(a==0 && b>0) { ref_haplotype[k] <- "pat"}
    if(a>0 && b==0) { ref_haplotype[k] <- "mat"}
    if(a==0 && b==length(thisSet_pat)) { ref_haplotype_strict[k] <- "pat"}
    if(a==length(thisSet_mat) && b==0) { ref_haplotype_strict[k] <- "mat"}

  }

  vg_RAF <- rep(NA, nrow(this_vg_bowtie))
  bowtie_RAF <- rep(NA, nrow(this_vg_bowtie))
  vg_RAF_strict <- rep(NA, nrow(this_vg_bowtie))
  bowtie_RAF_strict <- rep(NA, nrow(this_vg_bowtie))
  for (k in 1:nrow(this_vg_bowtie)) {
    if(!is.na(ref_haplotype[k]) && ref_haplotype[k]=="mat") {
      vg_RAF[k] <- this_vg_bowtie[k,"matCount"] / (this_vg_bowtie[k,"matCount"] + this_vg_bowtie[k,"patCount"])
      bowtie_RAF[k] <- this_vg_bowtie[k,"def_matCount"] / (this_vg_bowtie[k,"def_matCount"] + this_vg_bowtie[k,"def_patCount"])
    }
    if(!is.na(ref_haplotype[k]) && ref_haplotype[k]=="pat") {
      vg_RAF[k] <- this_vg_bowtie[k,"patCount"] / (this_vg_bowtie[k,"matCount"] + this_vg_bowtie[k,"patCount"])
      bowtie_RAF[k] <- this_vg_bowtie[k,"def_patCount"] / (this_vg_bowtie[k,"def_matCount"] + this_vg_bowtie[k,"def_patCount"])
    }
    if(!is.na(ref_haplotype_strict[k]) && ref_haplotype_strict[k]=="mat") {
      vg_RAF_strict[k] <- this_vg_bowtie[k,"matCount"] / (this_vg_bowtie[k,"matCount"] + this_vg_bowtie[k,"patCount"])
      bowtie_RAF_strict[k] <- this_vg_bowtie[k,"def_matCount"] / (this_vg_bowtie[k,"def_matCount"] + this_vg_bowtie[k,"def_patCount"])
    }
    if(!is.na(ref_haplotype_strict[k]) && ref_haplotype_strict[k]=="pat") {
      vg_RAF_strict[k] <- this_vg_bowtie[k,"patCount"] / (this_vg_bowtie[k,"matCount"] + this_vg_bowtie[k,"patCount"])
      bowtie_RAF_strict[k] <- this_vg_bowtie[k,"def_patCount"] / (this_vg_bowtie[k,"def_matCount"] + this_vg_bowtie[k,"def_patCount"])
    }

  }

  finalMat <- cbind(this_vg_bowtie, vg_RAF, bowtie_RAF, vg_RAF_strict, bowtie_RAF_strict)
  return(finalMat)

}


################################################################################
#This function determines, for each haplotype, the smallest MAF present in
#the haplotype, and adds it to the end of the column.
################################################################################
get_MAF_info <- function(this_comparison_Matrix, vep_d1, vep_d2) {
  #this_comparison_Matrix <- comparison_Matrix
  this_comparison_Matrix$uniqueVariants <- paste(this_comparison_Matrix$matVars, this_comparison_Matrix$patVars, sep="_")
  uniqueVariants <- unique(this_comparison_Matrix[,c("chr", "start", "end", "matVars", "patVars", "uniqueVariants")])
  uniqueVariants_gr <- makeGRangesFromDataFrame(uniqueVariants)

  all_vep <- rbind(vep_d1, vep_d2)
  all_vep_locations <- as.data.frame(matrix(nrow=nrow(all_vep), ncol=2, data=unlist(strsplit(as.character(all_vep$Location), split=":", fixed=T)), byrow=T))
  colnames(all_vep_locations) <- c("chr", "start")
  withDash <- grep("-", all_vep_locations$start)
  all_vep_locations[-withDash, "start"] <- paste(all_vep_locations[-withDash, "start"], all_vep_locations[-withDash, "start"], sep="-")
  corrected <- as.data.frame(matrix(nrow=nrow(all_vep_locations), ncol=2, data=unlist(strsplit(all_vep_locations[,2], split="-")), byrow=T))
  all_vep_locations <- cbind(all_vep_locations$chr, corrected)
  colnames(all_vep_locations) <- c("chr", "start", "end")
  all_vep_locations_gr <- makeGRangesFromDataFrame(all_vep_locations)

  theIntersect_1 <- as.data.frame(findOverlaps(uniqueVariants_gr, all_vep_locations_gr))

  min_maf <- rep(NA, nrow(uniqueVariants))
  uniquePlaces <- unique(theIntersect_1[,1])
  for (k in uniquePlaces) {
    if(which(uniquePlaces==k)%%1000==0) { print(paste(which(uniquePlaces==k), length(uniquePlaces)))}
    theSet <- theIntersect_1[theIntersect_1[,1]==k,]

    this_miniVep <- all_vep[theSet[,2],]
    theChrom <- all_vep_locations[theSet[1,2],1]


    thisSet_mat <- strsplit(uniqueVariants[theSet[,1],"matVars"], split=";")[[1]]
    thisSet_pat <- strsplit(uniqueVariants[theSet[,1],"patVars"], split=";")[[1]]
    thisSet_mat <- thisSet_mat[grep("HET", thisSet_mat)]
    thisSet_pat <- thisSet_pat[grep("HET", thisSet_pat)]

    allPresentVars <- c()
    for (l in 1:length(thisSet_mat)) {
      thisOne_mat <- strsplit(thisSet_mat[l], split=",")[[1]]
      thisOne_pat <- strsplit(thisSet_pat[l], split=",")[[1]]
      variant_first <- paste(theChrom, "_", thisOne_mat[1], "_", thisOne_mat[2], "/", thisOne_pat[2], sep="")
      variant_second <- paste(theChrom, "_", thisOne_mat[1], "_", thisOne_pat[2], "/", thisOne_mat[2], sep="")
      allPresentVars <- c(allPresentVars, variant_first, variant_second)
    }

    this_miniVep <- this_miniVep[this_miniVep$Uploaded_variation%in%allPresentVars,]
    this_miniVep$gnomad3_AF <- as.numeric(this_miniVep$gnomad3_AF)
    this_miniVep <- this_miniVep[!is.na(this_miniVep$gnomad3_AF),]
    if(nrow(this_miniVep)>0) {
      this_maf <- min(this_miniVep$gnomad3_AF)
      min_maf[k] <- this_maf
    }

  }

  uniqueVariants$minMAF <- min_maf
  theMatch <- match(this_comparison_Matrix$uniqueVariants, uniqueVariants$uniqueVariants)
  this_comparison_Matrix$minMAF <- uniqueVariants[theMatch,"minMAF"]

  this_comparison_Matrix <- this_comparison_Matrix[,!colnames(this_comparison_Matrix)%in%c("uniqueVariants")]
  return(this_comparison_Matrix)

}



################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################

args <- commandArgs(trailingOnly=T)
outDir <- args[1]
vg_directory <- args[2]
bowtie2_directory <- args[3]
variants_rep1 <- args[4]
variants_rep2 <- args[5]
pval_table_1 <- args[6]
pval_table_2 <- args[7]

#outDir <- "/cluster/home/bmoyers/Figures/BrainTF_AlleleSpecificBinding/Paper_Figures_2023September15/"
#vg_directory <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/ComparingBryanSteph/Donor_vg_test"
#bowtie2_directory <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/ComparingBryanSteph/Default_wContigs"
#variants_rep1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta"
#variants_rep2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences.fasta"
#pval_table_1 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#pval_table_2 <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/HaplotypePvalsTable_binom_vg_summedAcrossTissues.txt"
#vep_d1_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0001_annotated.tsv.gz"
#vep_d2_saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0002_annotated.tsv.gz"


#/cluster/home/jlawlor/belle

thePileups <- list.files(bowtie2_directory, pattern="pileup", full.names=T)

bowtie2_results <- list.files(bowtie2_directory, pattern="variantReads_counts", full.names=T)
bowtie2_results <- bowtie2_results[-grep("Control", bowtie2_results)]

vg_results <- list.files(vg_directory, pattern="variantReads_counts", full.names=T)
vg_results <- vg_results[-grep("Control", vg_results)]

variants_rep1 <- parseVariants(variants_rep1)
variants_rep2 <- parseVariants(variants_rep2)

all_experiments <- find_all_experiment_names(bowtie2_results)
all_experiments <- unique(all_experiments)


################################################################################
#I would like to remove any case which is significant in the Input.
#While this was a calculation that was downstream of this comparison,
#we learned that these regions seem to have technical problems, so best to
#remove them entirely.
#
#Identify those regions which produced a significant input for donors.
################################################################################

pval_table_1 <- read.table(pval_table_1, header=T, sep="\t", stringsAsFactors=F)

significant_input_d1 <- pval_table_1[pval_table_1[,"INPUT"]<=0.05,]


pval_table_2 <- read.table(pval_table_2, header=T, sep="\t", stringsAsFactors=F)

significant_input_d2 <- pval_table_2[pval_table_2[,"INPUT"]<=0.05,]



################################################################################
#I initially was restricting to cases which were significant in at least one
#case. However, I will also want to look at information on non-significant
#cases. So I really need to consider everything available, so long as
#there is a minimum depth of 6 in at least one of the cases, to reduce
#computational load.
################################################################################

comparison_Matrix <- data.frame(chr=c(), start=c(), end=c(), matCount=c(), patCount=c(),
  totalDepth=c(), pval=c(), def_matCount=c(), def_patCount=c(),
  def_totalDepth=c(), default_pval=c(), matVars=c(), patVars=c(), vg_RAF=c(), bowtie_RAF=c(),
  vg_RAF_strict=c(), bowtie_RAF_strict=c(), Experiment=c())


#saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6.txt", sep="")
#comparison_Matrix <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)


for (i in 1:length(all_experiments)) {
  print(paste(i, length(all_experiments), sep="  "))
  this_vg_set <- vg_results[grep(all_experiments[i], vg_results)]
  this_bowtie2 <- bowtie2_results[grep(all_experiments[i], bowtie2_results)]

  print("Calculating Significance...")
  this_vg_bowtie <- getSignificances(this_vg_set, this_bowtie2)

  print("Adding Variant Info...")
  this_vg_bowtie <- get_variant_descriptions(this_vg_bowtie, all_experiments[i], variants_rep1, variants_rep2)

  this_vg_bowtie <- this_vg_bowtie[grep("HET", this_vg_bowtie[,"matVars"]),]

  print("Adding RAF Info...")
  this_vg_bowtie <- get_refAlt_info(this_vg_bowtie, all_experiments[i], thePileups)

  this_vg_bowtie <- cbind(this_vg_bowtie, rep(all_experiments[i], nrow(this_vg_bowtie)))
  colnames(this_vg_bowtie)[ncol(this_vg_bowtie)] <- "Experiment"

  print("Done!")
  comparison_Matrix <- rbind(comparison_Matrix, this_vg_bowtie)

  saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6.txt", sep="")
  write.table(comparison_Matrix, saveFile, col.names=T, row.names=F, sep="\t")

}



saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6.txt", sep="")
write.table(comparison_Matrix, saveFile, col.names=T, row.names=F, sep="\t")
#comparison_Matrix <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)

#saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6.txt", sep="")
#comparison_Matrix <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)


################################################################################
################################################################################
#We want to add MAF information to these variants.
################################################################################
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

saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0001_annotated_reduced.tsv"
write.table(vep_d1, saveFile)


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

saveFile <- "/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/James_vep_results/5397-JL-0002_annotated_reduced.tsv"
write.table(vep_d2, saveFile)



comparison_Matrix <- get_MAF_info(comparison_Matrix, vep_d1, vep_d2)
saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6.txt", sep="")
write.table(comparison_Matrix, saveFile, col.names=T, row.names=F, sep="\t")


################################################################################
################################################################################
#Below this, analyzing the results and making figures.
################################################################################
################################################################################



################################################################################
#Setting up some initial information.
################################################################################

comparison_Matrix$vg_negLog10Pval <- -log10(comparison_Matrix$pval + 0.000000000000000001)
comparison_Matrix$bowtie_negLog10Pval <- -log10(comparison_Matrix$default_pval + 0.000000000000000001)




saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6.txt", sep="")
write.table(comparison_Matrix, saveFile, col.names=T, row.names=F, sep="\t")
#comparison_Matrix <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)



################################################################################
#Remove significant input regions.
################################################################################
region_names <- paste(comparison_Matrix[,1], ":", comparison_Matrix[,2], "-", comparison_Matrix[,3], sep="")
toRemove <- which(region_names%in%c(rownames(significant_input_d1), rownames(significant_input_d2)))

comparison_Matrix <- comparison_Matrix[-toRemove,]


saveFile <- paste(outDir, "Comparison_Matrix_FairRegionComparison_minDepth6_noSigInput.txt", sep="")
write.table(comparison_Matrix, saveFile, col.names=T, row.names=F, sep="\t")
#comparison_Matrix <- read.table(saveFile, header=T, sep="\t", stringsAsFactors=F)


a <- comparison_Matrix[which(comparison_Matrix$default_pval<=0.05),]
a <- a[a$minMAF>0.05,]
a <- a[a$minMAF<0.95,]
dim(a)
#a <- a[!is.na(a$bowtie_RAF),]
#dim(a[a$bowtie_RAF>0.5,])
a <- comparison_Matrix[which(comparison_Matrix$pval<=0.05),]
a <- a[a$minMAF>0.05,]
a <- a[a$minMAF<0.95,]
dim(a)




graphDF_1 <- comparison_Matrix[,c("vg_RAF", "vg_RAF_strict", "vg_negLog10Pval", "minMAF")]
graphDF_1$method <- rep("vg", nrow(graphDF_1))
colnames(graphDF_1) <- c("RAF", "RAF_strict", "negLog10Pval", "minMAF", "method")
graphDF_2 <- comparison_Matrix[,c("bowtie_RAF", "bowtie_RAF_strict", "bowtie_negLog10Pval", "minMAF")]
graphDF_2$method <- rep("bowtie", nrow(graphDF_2))
colnames(graphDF_2) <- c("RAF", "RAF_strict", "negLog10Pval", "minMAF", "method")
graphDF <- rbind(graphDF_1, graphDF_2)
graphDF <- as.data.frame(graphDF)

graphDF$binned_pval <- rep("0 to 1.3 (ns)", nrow(graphDF))
graphDF[graphDF$negLog10Pval>=2,"binned_pval"] <- "2 to 3"
graphDF[graphDF$negLog10Pval>=3,"binned_pval"] <- "3 to 4"
graphDF[graphDF$negLog10Pval>=4,"binned_pval"] <- "4 to 5"
graphDF[graphDF$negLog10Pval>=5,"binned_pval"] <- "5 to 6"
graphDF[graphDF$negLog10Pval>=6,"binned_pval"] <- "6 to 10"
graphDF[graphDF$negLog10Pval>=10,"binned_pval"] <- "10 to 15"
graphDF[graphDF$negLog10Pval>=15,"binned_pval"] <- "15 to 20"
graphDF[graphDF$negLog10Pval>=20,"binned_pval"] <- "20+"




################################################################################
#Reference allele Frequency Changes.
################################################################################


mini_graphDF <- graphDF
mini_graphDF <- mini_graphDF[mini_graphDF$negLog10Pval>= -log10(0.05),]
mini_graphDF <- mini_graphDF[!is.na(mini_graphDF$RAF),]
mini_graphDF <- mini_graphDF[mini_graphDF$binned_pval!="0 to 1.3 (ns)",]

saveFile_5 <- paste(outDir, "Supplemental_Figure_1.pdf", sep="")
density_RAF <- ggplot(mini_graphDF, aes(color=method, fill=method, x=RAF)) + geom_density(alpha=0.5) +
  theme_minimal() + theme(axis.title=element_text(size=25), axis.text=element_text(size=20)) + theme(legend.text=element_text(size=20), legend.title=element_text(size=25), legend.position="bottom") +
  theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Reference Allele Frequency") + ylab("Density")
ggsave(saveFile_5)

mini_graphDF_vg <- mini_graphDF[mini_graphDF$method=="vg",]
mini_graphDF_vg <- mini_graphDF_vg[!is.na(mini_graphDF_vg$RAF),]
nrow(mini_graphDF_vg[mini_graphDF_vg$RAF<0.5,]) / nrow(mini_graphDF_vg)
nrow(mini_graphDF_vg[mini_graphDF_vg$RAF>0.5,]) / nrow(mini_graphDF_vg)

mini_graphDF_bowtie <- mini_graphDF[mini_graphDF$method=="bowtie",]
mini_graphDF_bowtie <- mini_graphDF_bowtie[!is.na(mini_graphDF_bowtie$RAF),]
nrow(mini_graphDF_bowtie[mini_graphDF_bowtie$RAF<0.5,]) / nrow(mini_graphDF_bowtie)
nrow(mini_graphDF_bowtie[mini_graphDF_bowtie$RAF>0.5,]) / nrow(mini_graphDF_bowtie)


mini_graphDF <- mini_graphDF[!is.na(mini_graphDF$minMAF),]
mini_graphDF <- mini_graphDF[mini_graphDF$minMAF>0.05,]
mini_graphDF <- mini_graphDF[mini_graphDF$minMAF<0.95,]

saveFile_5 <- paste(outDir, "Figure_1B.pdf", sep="")
density_RAF <- ggplot(mini_graphDF, aes(color=method, fill=method, x=RAF)) + geom_density(alpha=0.5) +
  theme_minimal() + theme(axis.title=element_text(size=25), axis.text=element_text(size=20)) + theme(legend.text=element_text(size=20), legend.title=element_text(size=25), legend.position="bottom") +
  theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Reference Allele Frequency") + ylab("Density")
ggsave(saveFile_5)

mini_graphDF_vg <- mini_graphDF[mini_graphDF$method=="vg",]
nrow(mini_graphDF_vg[mini_graphDF_vg$RAF<0.5,]) / nrow(mini_graphDF_vg)
nrow(mini_graphDF_vg[mini_graphDF_vg$RAF>0.5,]) / nrow(mini_graphDF_vg)

mini_graphDF_bowtie <- mini_graphDF[mini_graphDF$method=="bowtie",]
nrow(mini_graphDF_bowtie[mini_graphDF_bowtie$RAF<0.5,]) / nrow(mini_graphDF_bowtie)
nrow(mini_graphDF_bowtie[mini_graphDF_bowtie$RAF>0.5,]) / nrow(mini_graphDF_bowtie)



################################################################################
#Now for all variants, including non-significant.
################################################################################

mini_graphDF <- graphDF
mini_graphDF <- mini_graphDF[!is.na(mini_graphDF$RAF),]

saveFile_5 <- paste(outDir, "Supplemental_Figure_2.pdf", sep="")
density_RAF <- ggplot(mini_graphDF, aes(color=method, fill=method, x=RAF)) + geom_density(alpha=0.5) +
  theme_minimal() + theme(axis.title=element_text(size=25), axis.text=element_text(size=20)) + theme(legend.text=element_text(size=20), legend.title=element_text(size=25), legend.position="bottom") +
  theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Reference Allele Frequency") + ylab("Density")
ggsave(saveFile_5)



################################################################################
#Next, make a table of the fates of these different variants.
################################################################################


default_classifier <- rep("Sig_Ref", nrow(comparison_Matrix))
default_classifier[comparison_Matrix$bowtie_RAF<0.5] <- "Sig_Alt"
default_classifier[comparison_Matrix$default_pval>0.05] <- "NonSig"


vg_classifier <- rep("Sig_Ref", nrow(comparison_Matrix))
vg_classifier[comparison_Matrix$vg_RAF<0.5] <- "Sig_Alt"
vg_classifier[comparison_Matrix$pval>0.05] <- "NonSig"

comparison_Matrix$default_classifier <- default_classifier
comparison_Matrix$vg_classifier <- vg_classifier


a <- table(default_classifier, vg_classifier)

theValues <- c("Sig_Ref", "Sig_Alt", "NonSig")

alt_graphDF <- matrix(nrow=3, ncol=3, data=NA)
for (i in 1:length(theValues)) {
  vg_subset <- comparison_Matrix[comparison_Matrix$vg_classifier==theValues[i],]
  for (j in 1:length(theValues)) {
    bowtie_subset <- vg_subset[vg_subset$default_classifier==theValues[j],]
    alt_graphDF[i,j] <- nrow(bowtie_subset)
  }
}
rownames(alt_graphDF) <- theValues
colnames(alt_graphDF) <- theValues

alt_graphDF <- alt_graphDF[,c("Sig_Ref", "NonSig", "Sig_Alt")]
alt_graphDF <- alt_graphDF[c("Sig_Ref", "NonSig", "Sig_Alt"),]

col_fun = colorRamp2(c(0, 5000, 10225), c("white", "purple", "red"))

heatmap_labels <- list(c("Sig_Ref"), c("NonSig", "Alt_Ref"))

saveFile_6 <- paste(outDir, "Figure_1C.pdf", sep="")
pdf(saveFile_6)
Heatmap(alt_graphDF, col = col_fun, column_title = "bowtie2 classification", row_title = "vg classification",
  name = "foo", cell_fun = function(j, i, x, y, width, height, fill)
        {grid.text(sprintf("%.0f", alt_graphDF[i, j]), x, y, gp = gpar(fontsize = 30, fontface = "bold"))},
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_title_side = "bottom",
  row_names_gp = grid::gpar(fontsize = 20),
  column_names_gp = grid::gpar(fontsize = 20),
  column_title_gp = gpar(fontsize = 30),
  row_title_gp = gpar(fontsize = 30),
  show_heatmap_legend=FALSE,
  cluster_columns=FALSE,
  cluster_rows=FALSE
)
dev.off()


################################################################################
#Prove that changes are due to increased read depth.
################################################################################

theMismatch <- which(default_classifier!=vg_classifier)
b <- comparison_Matrix[theMismatch,]


summary(comparison_Matrix$totalDepth-comparison_Matrix$def_totalDepth)
summary(b$totalDepth-b$def_totalDepth)

b$depthImprovement <- b$totalDepth-b$def_totalDepth


saveFile_5 <- paste(outDir, "Supplemental_Figure_3.pdf", sep="")
density_RAF <- ggplot(b, aes(x=depthImprovement)) + geom_density(alpha=0.5, fill="purple") +
  theme_minimal() + theme(axis.title=element_text(size=25), axis.text=element_text(size=20)) + theme(legend.text=element_text(size=20), legend.title=element_text(size=25), legend.position="bottom") +
  theme(axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Increase in # of reads using vg") + ylab("Density") + xlim(-20,40)
ggsave(saveFile_5)


#
