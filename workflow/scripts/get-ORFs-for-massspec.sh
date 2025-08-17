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

#---------------------------------------------------------------------------------------------------#
# Make different versions of ORF files (all, clust, clust_contained) for mass-spec peptide analysis #
#---------------------------------------------------------------------------------------------------#
{
   # Set up output folder
   mkdir ${OUT_FOLDER}
   
   # Basic header reformatting
   cd ${OUT_FOLDER}
   perl -pe 's/\;.*//; s/ GeneID=/|/;'  ../${snakemake_input[reclocus_cds_aa]}  \
      >  ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.faa
   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.faa \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.tab
   awkt '{out=""} $2 ~ /^M/ {out="start"} $2 ~ /\*$/{out=out"stop"} {print $1,out}' ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.tab \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.tab.startstop

   # Run CD-HIT collapse relative to largest sequence
   cd-hit \
      -M 10000 -aL 1 -c 1 -d 0 \
      -i ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.faa \
      -o ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust.faa
   
   # Additionally CD-HIT collapse relative to smallest sequence (to get ORFs contained fully within another ORF)
   cd-hit \
      -M 10000 -aS 1 -c 1 -d 0 \
      -i ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust.faa \
      -o ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_contained.faa

   # Reformat headers for 'all' ORFs
   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.faa \
      | perl -pe 's/\|/\t/; s/\*$//;' \
      | awkt '{print ">cu|"$1"_0|"$1" "$2" [0-0]\n"$3}' \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full_header.faa

   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full.faa \
      | perl -pe 's/\|/\t/; s/\*$//;' \
      | awk '!( $2 in a){a[$2] = out++} {print ">cu|"$1"_0|"$1" gene"a[$2]" [0-0]\n"$3}' \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_full_header_generic.faa

   # Reformat headers for 'clust' ORFs
   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust.faa \
      | perl -pe 's/\|/\t/; s/\*$//;' \
      | awkt '{print ">cu|"$1"_0|"$1" "$2" [0-0]\n"$3}' \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_header.faa

   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust.faa \
      | perl -pe 's/\|/\t/; s/\*$//;' \
      | awk '!( $2 in a){a[$2] = out++} {print ">cu|"$1"_0|"$1" gene"a[$2]" [0-0]\n"$3}' \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_header_generic.faa

   # Reformat headers for 'clust_contained' ORFs
   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_contained.faa \
      | perl -pe 's/\|/\t/; s/\*$//;' \
      | awkt '{print ">cu|"$1"_0|"$1" "$2" [0-0]\n"$3}' \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_contained_header.faa

   fasta2tabbed.pl ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_contained.faa \
      | perl -pe 's/\|/\t/; s/\*$//;' \
      | awk '!( $2 in a){a[$2] = out++} {print ">cu|"$1"_0|"$1" gene"a[$2]" [0-0]\n"$3}' \
      > ${snakemake_params[prefix]}_reference_reclocus_CDS_aa_clust_contained_header_generic.faa

} 2>&1 | tee -a ${snakemake_log[out]}

