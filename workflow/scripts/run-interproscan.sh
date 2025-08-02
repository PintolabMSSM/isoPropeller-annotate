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

#------------------------------#
# Run interproscan distributed #
#------------------------------#
{
   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}
   
   # Run interproscan
   cd ${OUT_FOLDER}
   perl -pe 's/\*$//' ../${snakemake_input[reclocus_cds_aa]} > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_iprscan.fa
   interproscan-distributed.sh ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_iprscan.fa 750 ${snakemake_params[lsf_allocation]} 30 -appl Pfam,PANTHER,SMART -t p -ms 100

   # Create a check file to let snakemake know that we launched all jobs
   touch ../${snakemake_output[interpro_launched]}

} 2>&1 | tee -a ${snakemake_log[out]}
