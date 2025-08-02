#!/usr/bin/env bash

###############
# ENVIRONMENT #
###############

# Set paths from snakemake parameters
export PATH="${snakemake_params[path]}:$PATH"
export PERL5LIB="${snakemake_params[perl5lib]}:$PERL5LIB"

# Set up cDNA cupcake and python path
CWD=`pwd`
cd ${snakemake_params[cdna_cupcake_path]}
python setup.py build   >   "${CWD}/${snakemake_log[build]}" 2>&1
python setup.py install >>  "${CWD}/${snakemake_log[build]}" 2>&1
cd ${CWD}
export PYTHONPATH="${snakemake_params[cdna_cupcake_path]}/sequence:${snakemake_params[cdna_cupcake_path]}/annotation:$PYTHONPATH"

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

#-----------------------------------------------------------#
# Patch sqanti3 report script to allow for extra categories #
#-----------------------------------------------------------#
{
   echo -e "\n ## Patch sqanti3 reporting script ## \n"

   # We patch the lines where the plot categories are defined with a provided patch file
   cp "${snakemake_params[sq3_path]}/utilities/report_qc/SQANTI3_report.R" ${OUT_FOLDER}/SQANTI3_report.R
   patch ${OUT_FOLDER}/SQANTI3_report.R < ${snakemake_params[sq3_fusion_patch]}

} 2>&1 | tee -a ${snakemake_log[out]}


#---------------------------#
# Run sqanti3 report script #
#---------------------------#
{
   echo -e "\n ## Run sqanti3 reporting script ## \n"
   
   cut -f 1-48 ${snakemake_input[sq3_classification]} \
      | perl -pe 's/^#//' \
      | awkt '$6~/^fusion_novel/{$6="fusion_novel"} $6~/^fusion_known/{$6="fusion_known"} {print $0}'  \
      > ${FILE_PREFIX}_patched
   
   cut -f 1-48 ${snakemake_input[sq3_classification]} \
      | perl -pe 's/^#//' \
      | awkt '$6~/^fusion_novel/{$6="fusion_novel"} $6~/^fusion_known/{$6="fusion_known"} $5=="exons" || $5>1  {print $0}'  \
      > ${FILE_PREFIX}_patched_multiexon
   
   cut -f 1-48 ${snakemake_input[sq3_classification]} \
      | perl -pe 's/^#//' \
      | awkt '$6~/^fusion_novel/{$6="fusion_novel"} $6~/^fusion_known/{$6="fusion_known"} $5=="exons" || $5==1 {print $0}'  \
      > ${FILE_PREFIX}_patched_monoexon
   
   # Run full report
   Rscript ${OUT_FOLDER}/SQANTI3_report.R \
      ${FILE_PREFIX}_patched \
      ${snakemake_input[sq3_junctions]} \
      ${snakemake_input[sq3_params]} \
      "${snakemake_params[sq3_path]}/utilities/" \
      False \
      both \
      || : &
   
   # Run multiexon report
   Rscript ${OUT_FOLDER}/SQANTI3_report.R \
      ${FILE_PREFIX}_patched_multiexon \
      ${snakemake_input[sq3_junctions]} \
      ${snakemake_input[sq3_params]} \
      "${snakemake_params[sq3_path]}/utilities/" \
      False \
      both \
      || : &
   
   # Run monoexon report
   Rscript ${OUT_FOLDER}/SQANTI3_report.R \
      ${FILE_PREFIX}_patched_monoexon \
      ${snakemake_input[sq3_junctions]} \
      ${snakemake_input[sq3_params]} \
      "${snakemake_params[sq3_path]}/utilities/" \
      False \
      both \
      || : &

   # Wait for all the backgrounded sqanti processes to finish
   wait

   # Ugly hack, but we added "|| :" to the previous commands to ignore errors when the report generation fails
   # Then we touch the html and pdf output so the pipeline will always continue even if there is an issue
   # with the report generation, which is not uncommon since it seems to be patched together with duct tape.
   touch ${FILE_PREFIX}_patched_SQANTI3_report.html
   touch ${FILE_PREFIX}_patched_SQANTI3_report.pdf

} 2>&1 | tee -a ${snakemake_log[out]}


#---------#
# Cleanup #
#---------#
{
   echo -e "\n ## Remove temporary files ## \n"

   rm -f ${OUT_FOLDER}/SQANTI3_report.R \
         ${FILE_PREFIX}_patched \
         ${FILE_PREFIX}_patched_multiexon \
         ${FILE_PREFIX}_patched_monoexon

}
