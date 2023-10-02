#!/bin/bash
#meme_pipeline.sh

################################################################################
# The script below was used to produce meme outputs for use in downstream scripts,
#     The script below is tailored to our local compute cluster's specific
#     architecture, and will need to be rewritten for your own compute environment's
#     specific needs if you would like to reproduce its results.  In particular,
#     the various sourcings will likely need to be removed, as well as the log_msg
#     commands.  One will need to provide or define their own TEMP_DIR.
#
# Expected to be in the user's PATH:
#     meme 5.3.3
#     bedtools 2.28.0
#     python 2.7.15
#     R 3.6.1
#     Restrict_Peak_Size_summitCenter.R , provided in code.
#     Get_Flanking_Peak_Regions.R , provided in code
#     Motif_Control_Tests.R , provided in code.
#     MEME_Motif_Header.txt , provided in data, for construction of final files.
#     Scripts for centering peaks and extracting flanking sequence,
#     The nullseq_generate.py script from the LSGKM suite.
#
# Required Environment Variables:
#     BED_FILE            - gzipped Bed file containing the appropriate IDR peaks for analysis.
#                             this will be restricted to the top 10k peaks (maximum) for
#                             speed and computational efficiency.
#     BFILE               - A background file containing relevant background percentages for
#                             each nucleotide. Default is mammal_homo_sapiens_1000_199.na.bfile
#     GENOME              - path to the hg38 genome fasta file.
#     MATSEQS_1, MATSEQS_2   - FASTA files containing all of the haplotype 1 variant sequences
#                              for donors 1 and 2, respectively
#     PATSEQS_1, PATSEQS_2   - FASTA files containing all of the haplotye 2 variant sequences
#                              for donors 1 and 2, respectively
#     NULLSEQ_INDICES     - path to the nullseq indices for hg38.
#
#
# Optional Environment Variables:
#     OUTPUT_DIR          - Change default output directory
#
#
################################################################################
#BED_FILE=/cluster/home/jloupe/Outputs_aquas/Apr2019_ChIPseq/ASH2L_CB_Output/peak/spp/idr/optimal_set/ASH2L_CB_Output_ppr.IDR0.05.filt.narrowPeak.gz
#BFILE="/gpfs/gpfs1/home/bmoyers/ENCODE_lsgkm_meme_SMRTS/Scripts/mammal_homo_sapiens_1000_199.na.bfile"
#GENOME="hg38"
#MATSEQS_1="/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences_maternal.fasta"
#PATSEQS_1="/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0001/AllVar_Sequences_paternal.fasta"
#MATSEQS_2="/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences_maternal.fasta"
#PATSEQS_2="/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/5397-JL-0002/AllVar_Sequences_paternal.fasta"
#OUTPUT_DIR=/cluster/home/bmoyers/BrainTF_AlleleSpecificBinding/dsvm_analyses/ASH2L_CB



if [ -z ${BFILE+x} ]; then
  BFILE=/cluster/home/bmoyers/ENCODE_lsgkm_meme_SMRTS/Scripts/mammal_homo_sapiens_1000_199.na.bfile
fi


### Export an appropriate python path for this analysis
#export PYTHONPATH=/gpfs/gpfs2/software/python2.7.2

### Load relevant modules
module load g/meme/5.3.3
module load g/bedtools/2.28.0
module load g/python/2.7.15
module load g/python/2.7-modules
module load g/R/3.6.1

TEMP_DIR=$OUTPUT_DIR/temp_meme
mkdir $TEMP_DIR

### Set up output files
BED_TEMP=$TEMP_DIR/Peaks_Temp.bed
BED_SORTED_TEMP1=$TEMP_DIR/Peaks_sorted_Temp1.bed
BED_TOP500_TEMP2=$TEMP_DIR/Peaks_top500_Temp2.bed
BED_SORTED_TEMP3=$TEMP_DIR/Peaks_Temp_removeProblems.bed




if file $BED_FILE | grep "gz"
  then zcat $BED_FILE > $BED_TEMP
  else cp $BED_FILE $BED_TEMP
fi



sort -r -n -k 7 $BED_TEMP | cut -f 1,2,3,7,10 > $BED_SORTED_TEMP1


grep -v "chrM" $BED_SORTED_TEMP1 | grep -v "random" | grep -v "EBV" | grep -v "chrUn" > $BED_SORTED_TEMP3
mv $BED_SORTED_TEMP3 $BED_SORTED_TEMP1


head -500 $BED_SORTED_TEMP1 > $BED_TOP500_TEMP2

Rscript Restrict_Peak_Size_summitCenter.R $BED_TOP500_TEMP2 50




NUM_PEAKS=$(< "$BED_SORTED_TEMP1" wc -l)

if [ $NUM_PEAKS -ge 1000 ]; then

  FASTA_TOP500=${TEMP_DIR}/Peaks_for_motif_calling.fasta
  MOTIF_OUT=${OUTPUT_DIR}/meme_results
  CENTRIMO_OUT=${OUTPUT_DIR}/centrimo_results

  BED_T1=${TEMP_DIR}/Peaks_T1.bed
  BED_C1=${TEMP_DIR}/Peaks_C1.bed
  BED_T2=${TEMP_DIR}/Peaks_T2.bed
  BED_C2=${TEMP_DIR}/Peaks_C2.bed


  FASTA_T1=${TEMP_DIR}/Peaks_T1.fasta
  FASTA_C1=${TEMP_DIR}/Peaks_C1.fasta
  FASTA_T2=${TEMP_DIR}/Peaks_T2.fasta
  FASTA_C2=${TEMP_DIR}/Peaks_C2.fasta


  FIMO_T1=${TEMP_DIR}/fimo_T1
  FIMO_C1=${TEMP_DIR}/fimo_C1
  FIMO_T2=${TEMP_DIR}/fimo_T2
  FIMO_C2=${TEMP_DIR}/fimo_C2


  fastaFromBed -fo $FASTA_TOP500 -fi ${GENOME} -bed $BED_TOP500_TEMP2

  head -1000 $BED_SORTED_TEMP1 | tail -500 > $BED_T1
  Rscript Restrict_Peak_Size_summitCenter.R $BED_T1 150
  fastaFromBed -fo $FASTA_T1 -fi ${GENOME} -bed $BED_T1

  python nullseq_generate.py -x 1 -m 1000 -r 1 -o $BED_C1 $BED_TOP500_TEMP2 hg38 ${NULLSEQ_INDICES}
  Rscript Restrict_Peak_Size.R $BED_C1 150
  fastaFromBed -fo $FASTA_C1 -fi ${GENOME} -bed $BED_C1

  tail -n +501 $BED_SORTED_TEMP1 > $BED_T2
  Rscript Restrict_Peak_Size_summitCenter.R $BED_T2 150
  fastaFromBed -fo $FASTA_T2 -fi ${GENOME} -bed $BED_T2

  Rscript Get_Flanking_Peak_Regions.R $BED_T2 $BED_C2 300
  fastaFromBed -fo $FASTA_C2 -fi ${GENOME} -bed $BED_C2



  meme -oc $MOTIF_OUT -bfile $BFILE -dna -mod zoops -nmotifs 5 -minw 6 -maxw 20 $FASTA_TOP500
  MEME_OUTPUT_FILE="${MOTIF_OUT}/meme.txt"


  centrimo --local -oc $CENTRIMO_OUT -bfile $BFILE $FASTA_T2 $MEME_OUTPUT_FILE
  CENTRIMO_OUTPUT_FILE="${CENTRIMO_OUT}/centrimo.tsv"


  fimo --oc $FIMO_T1 $MEME_OUTPUT_FILE $FASTA_T1
  fimo --oc $FIMO_C1 $MEME_OUTPUT_FILE $FASTA_C1
  fimo --oc $FIMO_T2 $MEME_OUTPUT_FILE $FASTA_T2
  fimo --oc $FIMO_C2 $MEME_OUTPUT_FILE $FASTA_C2

  PASSED_MOTIFS=${OUTPUT_DIR}/Passed_Motifs
  mkdir $PASSED_MOTIFS
  PARTIAL_MOTIFS=${OUTPUT_DIR}/Partial_Motifs
  mkdir $PARTIAL_MOTIFS
  FAILED_MOTIFS=${OUTPUT_DIR}/Failed_Motifs
  mkdir $FAILED_MOTIFS

  Rscript Motif_Control_Tests.R $MEME_OUTPUT_FILE $FIMO_T1/fimo.tsv $FIMO_C1/fimo.tsv $FIMO_T2/fimo.tsv $FIMO_C1/fimo.tsv $PASSED_MOTIFS $PARTIAL_MOTIFS $FAILED_MOTIFS $BED_FILE $BED_C1 $BED_T2 $BED_C2 $MOTIF_OUT ${CENTRIMO_OUTPUT_FILE}

  ALL_PASSED_MOTIFS=${OUTPUT_DIR}/All_Passed_Motifs.txt
  PARTIAL_PASSED_MOTIFS=${OUTPUT_DIR}/Partial_Passed_Motifs.txt
  ALL_MOTIFS=${OUTPUT_DIR}/All_Motifs.txt
  cat MEME_Motif_Header.txt $PASSED_MOTIFS/*.txt > $ALL_PASSED_MOTIFS
  cat MEME_Motif_Header.txt $PARTIAL_MOTIFS/*.txt > $PARTIAL_PASSED_MOTIFS
  cat MEME_Motif_Header.txt $PASSED_MOTIFS/*.txt $PARTIAL_MOTIFS/*.txt  $FAILED_MOTIFS/*.txt > $ALL_MOTIFS


  if  cat ${ALL_PASSED_MOTIFS} | grep "MOTIF" ; then

    fimo --oc ${OUTPUT_DIR}Donor1_Maternal_FIMO/ ${ALL_PASSED_MOTIFS} ${MATSEQS_1}

    fimo --oc ${OUTPUT_DIR}Donor1_Paternal_FIMO/ ${ALL_PASSED_MOTIFS} ${PATSEQS_1}

    fimo --oc ${OUTPUT_DIR}Donor2_Maternal_FIMO/ ${ALL_PASSED_MOTIFS} ${MATSEQS_2}

    fimo --oc ${OUTPUT_DIR}Donor2_Paternal_FIMO/ ${ALL_PASSED_MOTIFS} ${PATSEQS_2}


  fi

fi



rm -rf ${TEMP_DIR}

#Tests for motif passing were based on https://genome.cshlp.org/content/suppl/2012/08/22/22.9.1798.DC1/FigS1.pdf
