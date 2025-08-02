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

#--------------------------#
# Get isoform count matrix #
#--------------------------#
{
   echo -e "\n ## Get isoform count matrix ## \n"
   
   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}
   
   # Just copy the count file
   cat ${snakemake_input[niap_merge_counts]} > ${snakemake_output[iso_count_matrix]}

} 2>&1 | tee -a ${snakemake_log[out]}
