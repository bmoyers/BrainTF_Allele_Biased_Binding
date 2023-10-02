#Parsing_Sam_Reads_In_Variant_Regions.R

################################################################################
#This script is used to compile restrict sam files (which have been surjected from
#     a .gam mapping) to only those which overlap with our variant regions. This
#     script is called by DonorSpecificMapping.sh .
#
#Run under R 4.1.0, but can be used under other versions with GenomicRanges available.
#
#This program takes the following arguments:
# thisSam: The Path to the .sam file from which reads will be read.
# outputFile: The Path to the .sam file to which the reads in variant regions
#     should be written.
# allVarRegions: A FASTA file which contains all regions that contain a
#     non-reference base, whether heterozygous or homozygous. Produced by the
#     Building_fasta_variant_sequences.sh script.
#
################################################################################



################################################################################
################################################################################
#Load Libraries and define options.
################################################################################
################################################################################


R.Version()
.libPaths()

options(scipen="100")
library(GenomicRanges)



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
thisSam <- args[1]
outputFile <- args[2]
allVarRegions <- args[3]
#args <- c()


################################################################################
#Read in the variant regions and make a genomic ranges object.
################################################################################

allVarRegions <- readLines(allVarRegions)
allVarRegions <- allVarRegions[grep(">", allVarRegions)]

allMaternal <- allVarRegions[grep(":M:", allVarRegions)]
allPaternal <- allVarRegions[grep(":P:", allVarRegions)]
allVarBed <- gsub(">", "", allMaternal)
allVarBed <- matrix(nrow=length(allVarBed), ncol=2, data=unlist(strsplit(allVarBed, split=":M:")), byrow=T)
allVarBed <- matrix(nrow=nrow(allVarBed), ncol=3, data=unlist(strsplit(allVarBed[,1], split=":|-")), byrow=T)
allVarBed <- as.data.frame(allVarBed)
allVarBed[,2] <- as.numeric(as.character(allVarBed[,2]))
allVarBed[,3] <- as.numeric(as.character(allVarBed[,3]))

allMaternal <- matrix(nrow=length(allMaternal), ncol=2, data=unlist(strsplit(allMaternal, split=":M:")), byrow=T)
allPaternal <- matrix(nrow=length(allPaternal), ncol=2, data=unlist(strsplit(allPaternal, split=":P:")), byrow=T)
allVarBed[,4] <- allMaternal[,2]
allVarBed[,5] <- allPaternal[,2]
colnames(allVarBed) <- c("chr", "start", "end", "maternal", "paternal")

allVarGR <- makeGRangesFromDataFrame(allVarBed, keep.extra.columns=TRUE, ignore.strand=TRUE)
rm(allMaternal, allPaternal, allVarBed, allVarRegions)


################################################################################
#Read in the sam file and make a genomic ranges object.
################################################################################


allReads <- read.delim(thisSam, header=F, comment.char="@", sep="\t", stringsAsFactors=F)
allReads <- allReads[!is.na(allReads[,2]),]
allReads_chr <- as.character(allReads[,3])
allReads_start <- as.numeric(as.character(allReads[,4]))
allReads_nchar <- nchar(as.character(allReads[,10]))
allReads_stop <- as.numeric(as.character(allReads_start + allReads_nchar))
allReads_names <- as.character(allReads[,1])
allReads_GR <- data.frame(chr=allReads_chr, start=allReads_start, end=allReads_stop, name=allReads_names)
allReads_GR <- allReads_GR[!is.na(allReads_GR[,2]),]


allReads_GR <- makeGRangesFromDataFrame(allReads_GR, keep.extra.columns=TRUE, ignore.strand=TRUE)
rm(allReads_chr, allReads_start, allReads_nchar, allReads_stop, allReads_names)


################################################################################
#Identify any reads with overlaps.  Save the results.
################################################################################


theIntersect <- findOverlaps(allReads_GR, allVarGR)
theIntersect <- data.frame(theIntersect)
theIntersect <- unique(theIntersect[,1])
theIntersect <- allReads_GR[theIntersect]
theIntersect <- data.frame(theIntersect)
theIntersect <- theIntersect[,"name"]

rm(allReads_GR, allVarGR)

theIntersect <- allReads[allReads[,1]%in%theIntersect,]

write.table(theIntersect, outputFile, col.names=F, row.names=F, quote=F, sep="\t")
