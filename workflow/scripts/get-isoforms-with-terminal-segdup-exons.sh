#!/usr/bin/env bash

###############
# ENVIRONMENT #
###############

# Set paths from snakemake parameters
export PATH="${snakemake_params[path]}:$PATH"
export PERL5LIB="${snakemake_params[perl5lib]}:$PERL5LIB"

# Set aliases
shopt -s expand_aliases
alias awkt="awk -F '\t' -v OFS='\t'"
alias sortt="sort -t $'\t'"

# Copy snakemake variables to local variables for readability
OUT_FOLDER="${snakemake_params[output_dir]}"
FILE_PREFIX="${OUT_FOLDER}/${snakemake_params[prefix]}"


########
# MAIN #
########

#---------------------------------------------#
# Identify terminal exons that map in segdups #
#---------------------------------------------#
{
   echo -e "\n ## Prepare inputs ## \n"
   
   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}
   
   # Convert sqanti3 gtf file to bed format
   gtf2gff.pl -a transcript_id ${snakemake_input[sq3_gtf]}     > ${FILE_PREFIX}_corrected.tmp.gff
   gff2bed.pl ${FILE_PREFIX}_corrected.tmp.gff                 > ${FILE_PREFIX}_corrected.tmp.bed
   bed2intronexongff.pl -v 1 ${FILE_PREFIX}_corrected.tmp.bed  > ${FILE_PREFIX}_corrected.intronexon.gff
   
   # Collapse all reference genome transcripts based on any exon overlap on the same strand
   gtf-get-gene-regions.pl ${snakemake_params[refgenome_gtf]} > ${FILE_PREFIX}_reference-gene-regions.gtf
   
   # Overlap the gff file with segdups and gene regions and find isoforms with terminal exons in segdups
   for LEVEL in 1 2 3 4
   do
      find-mismapped-terminal-exons.pl \
         -i ${FILE_PREFIX}_corrected.intronexon.gff \
         -s ${snakemake_params[segdups_bed]} \
         -g ${FILE_PREFIX}_reference-gene-regions.gtf \
         -l ${LEVEL} \
         > ${FILE_PREFIX}_corrected_mismapping_level-${LEVEL}.txt
      
      awkt '($2=="no" && $3=="yes" && $6>0 && $7>100000) || ($2=="no" && $8=="yes" && $11>0 && $12>100000) {print $1}' \
         ${FILE_PREFIX}_corrected_mismapping_level-${LEVEL}.txt \
         > ${FILE_PREFIX}_level${LEVEL}-mismapped.txt
   done
   
   # Gather all mismappings
   cat ${FILE_PREFIX}_level?-mismapped.txt | sort | uniq > ${FILE_PREFIX}_mismapped-all.txt
   
   # Intersections to get collapsed bed file
   intersect-by-ids -ff ${FILE_PREFIX}_corrected.tmp.bed -fc 4 -if ${FILE_PREFIX}_mismapped-all.txt > ${FILE_PREFIX}_mismapped-all.bed
   collapseBed.pl ${FILE_PREFIX}_mismapped-all.bed > ${FILE_PREFIX}_mismapped-all_collapsed.bed
   awkt '{print $1,$2,$3,"mismaplocus_"count++" "$1":"$2"-"$3,$4}' ${FILE_PREFIX}_mismapped-all_collapsed.bed > ${FILE_PREFIX}_mismaplocus_regions.txt
   
   # Cleanup
   rm -f ${FILE_PREFIX}_corrected.tmp.gff \
         ${FILE_PREFIX}_corrected.tmp.bed \
         ${FILE_PREFIX}_corrected.intronexon.gff \
         ${FILE_PREFIX}_reference-gene-regions.gtf

} 2>&1 | tee -a ${snakemake_log[out]}

