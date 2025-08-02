#!/usr/bin/env bash

###############
# ENVIRONMENT #
###############

# Set paths from snakemake parameters
export PATH="${snakemake_params[path]}:${snakemake_params[path_maptools]}:$PATH"
export PERL5LIB="${snakemake_params[perl5lib]}:$PERL5LIB"

# Set aliases
shopt -s expand_aliases
alias awkt="awk -F '\t' -v OFS='\t'"
alias sortt="sort -t $'\t'"

# Copy snakemake variables to local variables for readability
OUT_FOLDER=${snakemake_params[output_dir]}
FILE_PREFIX="${OUT_FOLDER}/${snakemake_params[prefix]}"

########
# MAIN #
########

#-----------#
# Run cpat3 #
#-----------#
{
   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}
   
   # Convert sqanti3 gtf file to bed format
   gtf2gff.pl -a transcript_id ${snakemake_input[sq3_gtf]}  > ${FILE_PREFIX}_corrected.tmp.gff
   gff2bed.pl ${FILE_PREFIX}_corrected.tmp.gff              > ${FILE_PREFIX}_corrected.tmp.bed
   
   # Run cpatv3 (best-hit by probability)
   cpat.py \
      -d ${snakemake_params[cpat_prebuilt_logit]} \
      -r ${snakemake_params[refgenome_fasta]} \
      -x ${snakemake_params[cpat_prebuilt_hex]} \
      -g ${FILE_PREFIX}_corrected.tmp.bed \
      --top-orf=5 \
      --min-orf=${snakemake_params[minorf_cpat3]} \
      --best-orf=p \
      --log-file=${FILE_PREFIX}_corrected.cpatv3p18.log \
      -o ${FILE_PREFIX}_corrected.cpatv3p18 &
   
   # Run cpatv3 (best-hit by length)
   cpat.py \
      -d ${snakemake_params[cpat_prebuilt_logit]} \
      -r ${snakemake_params[refgenome_fasta]} \
      -x ${snakemake_params[cpat_prebuilt_hex]} \
      -g ${FILE_PREFIX}_corrected.tmp.bed \
      --top-orf=5 \
      --min-orf=${snakemake_params[minorf_cpat3]} \
      --best-orf=l \
      --log-file=${FILE_PREFIX}_corrected.cpatv3l18.log \
      -o ${FILE_PREFIX}_corrected.cpatv3l18 &

   # Wait for the cpat jobs to finish
   echo -e "\n ## Waiting for cpat jobs to finish ## \n"
   wait $(jobs -p)

   # Use the cpat annotation files to make a new GTF file based on CPAT best ORF annotations
   #echo -e "\n ## Cpat jobs finished, generate gtf files with CPAT ORFs ## \n"
   #CPAT_best2gtf.pl ${snakemake_input[sq3_gtf]} ${snakemake_output[cpat_prob_best]} ${snakemake_output[cpat_prob_best_gtf]}.tmp
   #CPAT_best2gtf.pl ${snakemake_input[sq3_gtf]} ${snakemake_output[cpat_leng_best]} ${snakemake_output[cpat_leng_best_gtf]}.tmp
   #awkt '$3 != "CDS"' ${snakemake_input[sq3_gtf]} >> ${snakemake_output[cpat_prob_best_gtf]}.tmp
   #awkt '$3 != "CDS"' ${snakemake_input[sq3_gtf]} >> ${snakemake_output[cpat_leng_best_gtf]}.tmp
   #sortt -k1,1 -k5,5n ${snakemake_output[cpat_prob_best_gtf]}.tmp > ${snakemake_output[cpat_prob_best_gtf]}.tmp.sorted
   #sortt -k1,1 -k5,5n ${snakemake_output[cpat_leng_best_gtf]}.tmp > ${snakemake_output[cpat_leng_best_gtf]}.tmp.sorted
   #gffread -T ${snakemake_output[cpat_prob_best_gtf]}.tmp.sorted > ${snakemake_output[cpat_prob_best_gtf]}
   #gffread -T ${snakemake_output[cpat_leng_best_gtf]}.tmp.sorted > ${snakemake_output[cpat_leng_best_gtf]}

   # Cleanup
   #echo -e "\n ## Do final cleanup ## \n"
   rm -f ${FILE_PREFIX}_corrected.tmp.gff \
         ${FILE_PREFIX}_corrected.tmp.bed #\
   #      ${snakemake_output[cpat_prob_best_gtf]}.tmp* \
   #      ${snakemake_output[cpat_leng_best_gtf]}.tmp*

} 2>&1 | tee -a ${snakemake_log[out]}

