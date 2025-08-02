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

#-------------#
# Prep folder #
#-------------#
{

   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir -p ${OUT_FOLDER}

} 2>&1 | tee -a ${snakemake_log[out]}


#-------------------------------#
# Identify nmd and poison exons #
#-------------------------------#
{
   
   # Run splice junction finder
   nmd_splice_junction_finder.pl \
      -i ${snakemake_input[reclocus_cds_gtf]} \
      -o ${FILE_PREFIX}_reference_reclocus_CDS \
      -s 1 \
      -m ${snakemake_params[nmdj_min_cds_len]} \
      -n ${snakemake_params[nmdj_distance]} \
      -t ${snakemake[threads]}

   # Run splice junction parser
   nmd_splice_junction_parser.pl \
      ${FILE_PREFIX}_reference_reclocus_CDS_nmd_sj_per_transcript.txt \
      ${FILE_PREFIX}_reference_reclocus_CDS_nmd_sj_per_transcript_parsed.txt 

} 2>&1 | tee -a ${snakemake_log[out]}

