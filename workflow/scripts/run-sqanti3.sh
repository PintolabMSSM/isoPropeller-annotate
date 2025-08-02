#!/usr/bin/env bash

###############
# ENVIRONMENT #
###############

# Set paths from snakemake parameters
export PATH="${snakemake_params[path]}:$PATH"
export PERL5LIB="${snakemake_params[perl5lib]}:$PERL5LIB"

# Set up cDNA cupcake and python path
#CWD=`pwd`
#cd ${snakemake_params[cdna_cupcake_path]}
#python setup.py build   >   "${CWD}/${snakemake_log[build]}" 2>&1
#python setup.py install >>  "${CWD}/${snakemake_log[build]}" 2>&1
#cd ${CWD}
export PYTHONPATH="${snakemake_params[cdna_cupcake_path]}/sequence:${snakemake_params[cdna_cupcake_path]}/annotation:$PYTHONPATH"

# Set aliases
shopt -s expand_aliases
alias awkt="awk -F '\t' -v OFS='\t'"
alias sortt="sort -t $'\t'"

# Copy snakemake variables to local variables for readability
FILE_PREFIX=${snakemake_params[prefix]}
OUT_FOLDER=${snakemake_params[output_dir]}

########
# MAIN #
########

#----------------#
# Prepare inputs #
#----------------#
{   
   # Run sqanti3
   rm -rf ${OUT_FOLDER}
   python ${snakemake_params[sqanti3_qc]} \
      -t ${snakemake[threads]} \
      -d ${OUT_FOLDER} \
      -o ${FILE_PREFIX} \
      --min_ref_len 1 \
      --report pdf \
      --coverage         ${snakemake_params[intron_coverage]} \
      --CAGE_peak        ${snakemake_params[cage_peaks_allmerge]} \
      --polyA_peak       ${snakemake_params[polya_bed]} \
      --polyA_motif_list ${snakemake_params[polya_motifs]} \
      ${snakemake_input[isocollapse_gtf]} \
      ${snakemake_params[refgenome_gtf]} \
      ${snakemake_params[refgenome_fasta]} \
      || :

} 2>&1 | tee -a ${snakemake_log[out]}
