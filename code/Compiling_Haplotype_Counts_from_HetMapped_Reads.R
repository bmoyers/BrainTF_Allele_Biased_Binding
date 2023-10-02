#!/usr/bin/env R
#Compiling_Haplotype_Counts_from_HetMapped_Reads.R


################################################################################
#The following script determines the number of reads mapping to maternal and paternal
#     haplotypes of a given genome in a given experiment.  It is assumed that there are
#     multiple potential mappings for each read, some of which may be aligned equally well.
#     The script parses, for each read, what the maximum alignment score is and which mappings
#     produce the maximum alignment score.  Based on where these equally-good mappings
#     occur, it records either a full read mapping to the appropriate location, a half-read
#     mapping to each of the maternal and paternal haplotype, or throws the read out.  The
#     output is two tables which record the total number of reads mapped to each heterozygous
#     location in each of the paternal and maternal cases.
#
#Run under R 3.4.3, but can be used under any other version that has stringdist,
#     genomicRanges, and BioStrings packages available.
#
#This program takes the following arguments:
# hetSam: a file in the sam format which includes only those reads which map
#     to a heterozygous in the donor-specific genome.  This should
#     be automatically produced as part of the script DonorSpecificMapping.sh
# allVarSeqs: A FASTA file which contains all regions that contain a
#     non-reference base, whether heterozygous or homozygous. Produced by the
#     Building_fasta_variant_sequences.sh script.
# outputBase: A baseline path to which the outputs should be saved.
#     Will have "maternal.txt" and "paternal.txt" appended for each of the
#     two files. Note that these names are arbitrary, and simply refer to
#     two haplotypes.  We have not confirmed parent of origin.
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
#Define Functions.
################################################################################
################################################################################


#The following function creates substrings of the offered sequence and vector
kmerSubstrings <- function(theVec, theSeq) {
  return(substr(theSeq, theVec[1], theVec[2]))
}

################################################################################
################################################################################
#Begin Script.
################################################################################
################################################################################

R.Version()

library(stringdist)
library(GenomicRanges)
library(Biostrings)
options("scipen"=100)

args <- commandArgs(trailingOnly=T)
hetSam <- args[1]
allVarSeqs <- args[2]
outputBase <- args[3]



################################################################################
#Set up a matrix and genomic ranges for the variant regions.
################################################################################

allVarSeqs <- readLines(allVarSeqs)
allVarTitles <- allVarSeqs[grep(">", allVarSeqs)]

matVarTitles <- allVarTitles[grep(":M:", allVarTitles)]
patVarTitles <- allVarTitles[grep(":P:", allVarTitles)]
matSeqs <- which(allVarSeqs%in%matVarTitles)
patSeqs <- which(allVarSeqs%in%patVarTitles)

matVarTitles <- matrix(nrow=length(matVarTitles), ncol=2, data=unlist(strsplit(matVarTitles, split=":M:")), byrow=T)
patVarTitles <- matrix(nrow=length(patVarTitles), ncol=2, data=unlist(strsplit(patVarTitles, split=":P:")), byrow=T)

matVarTitles <- cbind(matVarTitles, allVarSeqs[matSeqs+1])
patVarTitles <- cbind(patVarTitles, allVarSeqs[patSeqs+1])



locMat <- gsub(">", "", matVarTitles[,1])
locMat <- gsub(":M:", "", locMat)
locMat <- matrix(nrow=length(locMat), ncol=3, data=unlist(strsplit(locMat, split="-|:")), byrow=T)

seqMat <- data.frame(chr=locMat[,1], start=locMat[,2], end=locMat[,3], maternal=matVarTitles[,3], paternal=patVarTitles[,3])
seqMat[,1] <- as.character(seqMat[,1])
seqMat[,2] <- as.numeric(as.character(seqMat[,2]))
seqMat[,3] <- as.numeric(as.character(seqMat[,3]))
seqMat[,4] <- as.character(seqMat[,4])
seqMat[,5] <- as.character(seqMat[,5])

seqMat <- seqMat[order(seqMat[,1], seqMat[,2], seqMat[,3]),]
seqMat <- unique(seqMat)

rm(allVarSeqs, allVarTitles, locMat, matSeqs, matVarTitles, patSeqs, patVarTitles)

################################################################################
#Read in the sam file and create a genomic ranges for it.
################################################################################

totalRows <- system(paste("wc -l ", hetSam, sep=""), intern=T)
totalRows <- as.numeric(strsplit(totalRows, split=" ")[[1]][1])
theSam <- read.delim(hetSam, header=F, sep="\t", stringsAsFactors=F, nrow=totalRows)
samMat <- data.frame(chr=as.character(theSam[,3]), start=as.numeric(as.character(theSam[,4])), end=rep(NA, nrow(theSam)), readSeq=as.character(theSam[,10]))
samMat[,"end"] <- samMat[,"start"] + nchar(as.character(samMat[,"readSeq"]))
samMat[,1] <- as.character(samMat[,1])
samMat[,2] <- as.numeric(as.character(samMat[,2]))
samMat[,3] <- as.numeric(as.character(samMat[,3]))
samMat[,4] <- as.character(samMat[,4])

samMat <- samMat[order(samMat[,1], samMat[,2]),]

rm(theSam)

seqMat_complete <- seqMat
#samMat_safe <- samMat

samMat_gr <- makeGRangesFromDataFrame(samMat, ignore.strand=T)
seqMat_gr <- makeGRangesFromDataFrame(seqMat, ignore.strand=T)




##########################################################
#Get an initial intersect of this samFile with the haplotypes.
##########################################################
theIntersect <- as.data.frame(findOverlaps(samMat_gr, seqMat_gr, type="within"))



##########################################################
#Next, identify all of the hetSam files for the same TF,
#but in any tissue.
##########################################################

dir_hetSam <- strsplit(hetSam, split="/")[[1]]
relevant_TF <- dir_hetSam[length(dir_hetSam)]
relevant_TF <- strsplit(relevant_TF, split="_")[[1]][2]
dir_hetSam <- dir_hetSam[1:(length(dir_hetSam)-1)]
dir_hetSam <- paste(dir_hetSam, collapse="/")

all_hetSams <- list.files(dir_hetSam, pattern=relevant_TF)
all_hetSams <- all_hetSams[grep("hetsOnly", all_hetSams)]
all_hetSams <- paste(dir_hetSam, all_hetSams, sep="/")
all_hetSams <- all_hetSams[all_hetSams!=hetSam]

allIntersects <- theIntersect

for (i in 1:length(all_hetSams)) {
  print(paste(i, length(all_hetSams)))
  this_totalRows <- system(paste("wc -l ", all_hetSams[i], sep=""), intern=T)
  this_totalRows <- as.numeric(strsplit(this_totalRows, split=" ")[[1]][1])

  this_theSam <- read.delim(all_hetSams[i], header=F, sep="\t", stringsAsFactors=F, nrow=this_totalRows)
  this_samMat <- data.frame(chr=as.character(this_theSam[,3]), start=as.numeric(as.character(this_theSam[,4])), end=rep(NA, nrow(this_theSam)), readSeq=as.character(this_theSam[,10]))
  this_samMat[,"end"] <- this_samMat[,"start"] + nchar(as.character(this_samMat[,"readSeq"]))
  this_samMat[,1] <- as.character(this_samMat[,1])
  this_samMat[,2] <- as.numeric(as.character(this_samMat[,2]))
  this_samMat[,3] <- as.numeric(as.character(this_samMat[,3]))
  this_samMat[,4] <- as.character(this_samMat[,4])

  this_samMat <- this_samMat[order(this_samMat[,1], this_samMat[,2]),]

  this_samMat_gr <- makeGRangesFromDataFrame(this_samMat, ignore.strand=T)

  this_theIntersect <- as.data.frame(findOverlaps(this_samMat_gr, seqMat_gr, type="within"))

  allIntersects <- rbind(allIntersects, this_theIntersect)

  rm(this_totalRows, this_theSam, this_samMat, this_samMat_gr, this_theIntersect)

}



##########################################################
#Remove haplotypes that will be uninformative.  This is based on 2 factors.
#
#First, check to see which haplotypes have fewer than 6 reads in this
#experiment of interest.  Taht constitutes our candidate set of "don't even bother"
#
#Second, check to see which haplotypes have at least N reads, where N is the
#number of reads necessary to achieve significance (after BH p-val correction)
#in a sum of experiment, where the n for the correction is equal to the
#total number of haplotype regions.
##########################################################

#First, determine the minimum number of reads necessary.
readNumVector <- 1:100
pvalVector <- rep(1, 100)
for (i in 1:length(readNumVector)) {
  pvalVector[i] <- p.adjust(binom.test(0,readNumVector[i],alt="two.sided")$p.value,nrow(seqMat), method="BH")
}
minReads <- which(pvalVector<0.05)
minReads <- min(minReads)



intTable_1 <- table(theIntersect[,2])
candidates_toRemove <- as.numeric(names(intTable_1[intTable_1<=5]))

intTable_2 <- table(allIntersects[,2])
def_toKeep <- as.numeric(names(intTable_2[intTable_2>=minReads]))

toRemove <- candidates_toRemove[!candidates_toRemove%in%def_toKeep]



seqMat <- seqMat[-toRemove,]
seqMat_gr <- makeGRangesFromDataFrame(seqMat, ignore.strand=T)


newIntersect <- as.data.frame(findOverlaps(samMat_gr, seqMat_gr, type="within"))

seqMat <- seqMat[unique(newIntersect[,2]),]
samMat <- samMat[unique(newIntersect[,1]),]


rm(all_hetSams, allIntersects, candidates_toRemove, def_toKeep, dir_hetSam, intTable_1, intTable_2, minReads, newIntersect, pvalVector, readNumVector, relevant_TF, samMat_gr, seqMat_gr, theIntersect, toRemove, totalRows)



seqMat_gr <- makeGRangesFromDataFrame(seqMat, ignore.strand=T)
samMat_gr <- makeGRangesFromDataFrame(samMat, ignore.strand=T)
newIntersect <- as.data.frame(findOverlaps(samMat_gr, seqMat_gr, type="within"))



##########################################################
#Determine overall sequence similarity between alternate mapping possibilities,
#and assign reads as appropriate.
##########################################################


miniTitles <- paste(seqMat[,1], ":", seqMat[,2], "-", seqMat[,3], sep="")
miniPatCounts <- matrix(nrow=length(miniTitles), ncol=1, data=0)
miniMatCounts <- matrix(nrow=length(miniTitles), ncol=1, data=0)

theCol <- strsplit(outputBase, split="/")[[1]]
theCol <- strsplit(theCol[length(theCol)], split=".filt.nodup.srt")[[1]][1]
rownames(miniMatCounts) <- miniTitles
rownames(miniPatCounts) <- miniTitles
colnames(miniMatCounts) <- theCol
colnames(miniPatCounts) <- theCol


time1 <- Sys.time()
for (i in 1:nrow(seqMat)) {
#for (i in 1:1000) {
  if(i%%1000==0) {
    time2 <- Sys.time()
    print(paste(i, nrow(seqMat), difftime(time2, time1, units="mins"), sep="  "))
    time1 <- time2
  }
  t1 <- as.character(seqMat[i,4])
  t2 <- as.character(seqMat[i,5])

  total_mat <- 0
  total_pat <- 0

  miniSamMat <- samMat[newIntersect[newIntersect[,2]==i,1],]
  miniSamMat <- miniSamMat[order(nchar(miniSamMat[,"readSeq"])),]
  uniqLens <- unique(nchar(miniSamMat[,"readSeq"]))
  for (j in 1:length(uniqLens)) {
    this_miniSamMat <- miniSamMat[nchar(miniSamMat[,"readSeq"])==uniqLens[j],]

    theArray1 <- cbind(1:(nchar(t1)-(uniqLens[j]-1)), uniqLens[j]:nchar(t1))
    allSubstrings1 <- apply(X=theArray1, MARGIN=1, FUN=kmerSubstrings, theSeq=t1)

    theArray2 <- cbind(1:(nchar(t2)-(uniqLens[j]-1)), uniqLens[j]:nchar(t2))
    allSubstrings2 <- apply(X=theArray2, MARGIN=1, FUN=kmerSubstrings, theSeq=t2)

    for (k in 1:nrow(this_miniSamMat)) {
      q <- this_miniSamMat[k,"readSeq"]
      q_rev <- as.character(reverseComplement(DNAString(q)))

      sim1_fwd <- max(stringsim(q, allSubstrings1))
      sim1_rev <- max(stringsim(q_rev, allSubstrings1))
      sim2_fwd <- max(stringsim(q, allSubstrings2))
      sim2_rev <- max(stringsim(q_rev, allSubstrings2))

      sim1 <- max(sim1_fwd, sim1_rev)
      sim2 <- max(sim2_fwd, sim2_rev)
      if(sim1==sim2) {
        total_mat <- total_mat+0.5
        total_pat <- total_pat+0.5
      }
      if(sim1>sim2) {
        total_mat <- total_mat+1
      }
      if(sim1<sim2) {
        total_pat <- total_pat+1
      }

    }

  }
  miniMatCounts[i,1] <- total_mat
  miniPatCounts[i,1] <- total_pat

}


mini_counts <- cbind(miniMatCounts, miniPatCounts)






##########################################################
#Write the tables compiling depth for each haplotype.
##########################################################


theTitles <- paste(seqMat_complete[,1], ":", seqMat_complete[,2], "-", seqMat_complete[,3], sep="")
patCounts <- matrix(nrow=length(theTitles), ncol=1, data=0)
matCounts <- matrix(nrow=length(theTitles), ncol=1, data=0)

theCol <- strsplit(outputBase, split="/")[[1]]
theCol <- strsplit(theCol[length(theCol)], split=".filt.nodup.srt")[[1]][1]
rownames(matCounts) <- theTitles
rownames(patCounts) <- theTitles
colnames(matCounts) <- theCol
colnames(patCounts) <- theCol

#Need to incorporate the mini[M|P]atCounts into these tables...

theMatch <- match(rownames(mini_counts), rownames(matCounts))

matCounts[theMatch,] <- mini_counts[,1]
patCounts[theMatch,] <- mini_counts[,2]




output_maternal <- paste(outputBase, "_maternal.txt", sep="")
output_paternal <- paste(outputBase, "_paternal.txt", sep="")

write.table(matCounts, output_maternal, col.names=T, row.names=T, sep="\t", quote=F)
write.table(patCounts, output_paternal, col.names=T, row.names=T, sep="\t", quote=F)


##########################################################
#Calculate raw p-values of imbalance, save a table of these.
##########################################################


thePvals <- cbind(rep(1, nrow(matCounts)))
colnames(thePvals) <- gsub("_variantReads_counts", "", colnames(matCounts))
rownames(thePvals) <- rownames(matCounts)

theDepths <- cbind(rowSums(cbind(patCounts[,1], matCounts[,1])))
test_vector <- floor(apply(cbind(patCounts[,1], matCounts[,1]), 1, max))
to_replace <- which(theDepths[,1]>1)
bt <- function(a, b, p = 0.5) { format(binom.test(a, b, 0.5, alternative="two.sided")$p.value, scientific=F) }
thePvals_replace <- mapply(bt, test_vector[to_replace], theDepths[to_replace,1])
thePvals[to_replace,] <- thePvals_replace

thePvals <- as.data.frame(thePvals)
thePvals[,1] <- as.numeric(as.character(thePvals[,1]))

finalDestination <- paste(gsub("_counts", "", outputBase), "_binomial_pval.txt", sep="")

write.table(thePvals, finalDestination, row.names=T, col.names=T, sep="\t", quote=F)



#bsub -n 3 -R rusage[mem=24000] -We 100:00 -q c7normal -o /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/ComparingBryanSteph/TempTest/ALT5_CTCF_AMY_Oct2018_1224_bsubLog_recount.bsubLog.txt -J recounts_vg "module load g/R/3.4.3; Rscript /gpfs/gpfs1/home/bmoyers/Scripts/BrainTF_AlleleSpecificBinding/CompilingHaplotypeReads_multipleMappings_vg_v10_2020July24.R /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/MappedReads_vg/CTCF_AMY_Oct2018_1224.filt.nodup.srt.hetsOnly.sam /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences.fasta /gpfs/gpfs1/home/bmoyers/BrainTF_AlleleSpecificBinding/ComparingBryanSteph/TempTest/ALT5_CTCF_AMY_Oct2018_1224_variantReads_counts"
