#!/usr/bin/env bash

###############
# ENVIRONMENT #
###############

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

#-------------------------------------------------------------------#
# Run PoGo with a provided peptide map to check for new ORF support #
#-------------------------------------------------------------------#
{
   # NOTE: PoGo must be built with gcc/4.9.3 or greater to enable proper regexp support.
   
   # Set up output folder
   mkdir ${OUT_FOLDER}
   
   # Create a PoGo compatible fasta file
   gtf-count-attributes.pl -a transcript_id,gene_id ${snakemake_input[gtf_stopfix]} \
      | cut -f 1,2 \
      > ${FILE_PREFIX}_tid-gid.mapping
   perl -pe 's/ .*//' ${snakemake_input[reclocus_cds_aa]} \
      | fasta2tabbed.pl - \
      | perl -pe 's/\*$//' \
      > ${FILE_PREFIX}_ORF.tab
   join-by-ids ${FILE_PREFIX}_tid-gid.mapping ${FILE_PREFIX}_ORF.tab \
      | awkt '$1 !~ /^#/ {print ">transcript:" $1 " gene:" $2 "\n" $3}' \
      > ${FILE_PREFIX}_PoGo_input.fasta
   rm -f ${FILE_PREFIX}_tid-gid.mapping ${FILE_PREFIX}_ORF.tab
   
   # Get PoGo compatible gtf file
   cp ${snakemake_input[gtf_stopfix]} ${FILE_PREFIX}_PoGo_input.gtf
   
   # Copy the PoGo peptide list
   cp ${snakemake_params[pogo_peptides]} ${FILE_PREFIX}_PoGo_mm0.txt
   cp ${snakemake_params[pogo_peptides]} ${FILE_PREFIX}_PoGo_mm1.txt
   cp ${snakemake_params[pogo_peptides]} ${FILE_PREFIX}_PoGo_mm2.txt
   
   # Run PoGo at two different mismatch (mm) thresholds
   for i in 0 1
   do
      ${snakemake_params[pogo_bin]} \
         -fasta ${FILE_PREFIX}_PoGo_input.fasta \
         -gtf   ${FILE_PREFIX}_PoGo_input.gtf \
         -in    ${FILE_PREFIX}_PoGo_mm${i}.txt \
         -idv 1 \
         -mm ${i}
   done

} 2>&1 | tee -a ${snakemake_log[out]}

