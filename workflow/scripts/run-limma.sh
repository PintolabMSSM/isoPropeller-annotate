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

########
# MAIN #
########

#-------------#
# Prep folder #
#-------------#
{

   # Set up output folder
   mkdir -p ${OUT_FOLDER}

   # Strip out decimals from the raw count matrix to conform to expected inputs
   perl -pe 's/(\t\d+)\.\d+/$1/g; s/#TranscriptID\t//;' \
      ${snakemake_input[iso_count_matrix]} \
      > ${OUT_FOLDER}/exp-raw-counts.txt

   # Gather the sample IDs from the count file
   numbered-header ${OUT_FOLDER}/exp-raw-counts.txt \
      | awk '{print $1"\t"$1}' \
      | perl -pe 's/\tflcount_/\t/' \
      > ${OUT_FOLDER}/exp-design.ids

   # Match the sample IDs back to the master design file specified in config.yaml
   sorted-intersect-by-ids \
      -ff ${snakemake_params[exp_design_file]} \
      -if ${OUT_FOLDER}/exp-design.ids \
      | perl -pe 's/^#//' \
      > ${OUT_FOLDER}/exp-design.txt
   
   # Convert the NIAP classification file into an anno file
   cut-by-ids -f ${snakemake_input[reclocus_cds_trns_txt]} \
      transcript_id chr length gene_name gene_type  \
      | grep -vP '^transcript_id' \
      | awkt '$4==""{$4="NovelGene"} $5==""{$5="unknown_type"} {print $0}' \
      > ${OUT_FOLDER}/exp.anno
   add-header ${OUT_FOLDER}/exp.anno GeneID Chr Length Symbol type_of_gene

   # Copy the contrast file that goes with the design file
   cp ${snakemake_params[exp_contrasts_file]} ${OUT_FOLDER}/exp-contrasts.txt

   # Copy the models file that goes with the design file
   cp ${snakemake_params[exp_models_file]} ${OUT_FOLDER}/exp-models.txt
   
   # Prepare tissue-specific gene list to produce the highlighted genes plot
   cut -f 2,3 ${snakemake_params[exp_highlight_file]} > ${OUT_FOLDER}/exp-highlight.tmp
   join-by-ids -1 4 ${OUT_FOLDER}/exp.anno ${OUT_FOLDER}/exp-highlight.tmp | cut -f 1,6 | sort | uniq > ${OUT_FOLDER}/exp-highlight.txt

   # Cleanup
   rm -f ${OUT_FOLDER}/exp-design.ids ${OUT_FOLDER}/exp-highlight.tmp

} 2>&1 | tee -a ${snakemake_log[out]}


#-----------#
# Run limma #
#-----------#
{
   # Change into output dir
   cd ${OUT_FOLDER}
   
   # Run limma
   for PVAL in 0.0005
   do    
      run-limma \
         -i exp-raw-counts.txt \
         -d exp-design.txt \
         -c exp-contrasts.txt \
         -m exp-models.txt \
         -a exp.anno \
         -e exp-highlight.txt \
         -j protein_coding,lncRNA,processed_pseudogene,unprocessed_pseudogene,misc_RNA,transcribed_unprocessed_pseudogene,transcribed_processed_pseudogene \
         -s human \
         -p ${PVAL} \
         -F 0.5,10 \
         -x DisplayID 
   done

   # Create the checkpoint file
   touch ${snakemake_params[prefix]}_limma.checkpoint
   
} 2>&1 | tee -a ${snakemake_log[out]}

