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
OUT_FOLDER=${snakemake_params[output_dir]}
FILE_PREFIX="${OUT_FOLDER}/${snakemake_params[prefix]}"

########
# MAIN #
########

#------------------------#
# Get a splicechain file #
#------------------------#
{
   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}
   
   # Convert the patched sqanti bedfile to 
   bed2splicechains.pl ${snakemake_input[bed_classpatched]} \
      | sort -t : -k1,1 -k2,2n \
      > ${FILE_PREFIX}_patched_splicechains.txt

} 2>&1 | tee -a ${snakemake_log[out]}

