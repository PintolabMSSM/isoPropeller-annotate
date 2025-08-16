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

#-------------#
# Prep folder #
#-------------#

{

   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}

} 2>&1 | tee -a ${snakemake_log[out]}


#------------------------#
# Run isoPropeller merge #
#------------------------#
{
   
   # Annotation
   isoPropeller_annotate -q ${snakemake_input[isocollapse_gtf]} -r ${snakemake_params[refgenome_isop_gtf]} -o ${snakemake_output[ref_gtf]} -t ${snakemake[threads]} -e ${snakemake_input[tss_bed]}
   isoPropeller_annotate -q ${snakemake_input[isocollapse_gtf]} -r ${snakemake_params[refgenome_isop_gtf]} -o ${snakemake_output[ref_im_gtf]} -t ${snakemake[threads]} -m -e ${snakemake_input[tss_bed]}

   # Tabulation   
   for attribute in `echo gene_name gene_type status asm_gene_id ref_transcript_id ref_gene_id ref_gene_name ref_gene_type`; do echo $attribute; done > ${FILE_PREFIX}_temp_attribute_list.txt
   gtf2summary.pl -i ${snakemake_output[ref_gtf]} -o ${FILE_PREFIX}_reference -a ${FILE_PREFIX}_temp_attribute_list.txt -g ${snakemake_params[refgenome_fasta]} -j ${snakemake_params[intron_coverage]} -r -t ${snakemake[threads]}
   gtf2summary.pl -i ${snakemake_output[ref_im_gtf]} -o ${FILE_PREFIX}_reference_im -a ${FILE_PREFIX}_temp_attribute_list.txt -g ${snakemake_params[refgenome_fasta]} -j ${snakemake_params[intron_coverage]} -r -t ${snakemake[threads]}

   # Locus reconstruction
   sed 's/#TranscriptID/transcript_id/' ${snakemake_input[isop_merge_counts]} > ${FILE_PREFIX}_temp_count.txt
   merge2tables.pl -t1 ${snakemake_output[ref_trns_txt]} -c1 0 -t2 ${FILE_PREFIX}_temp_count.txt -c2 0 -o ${FILE_PREFIX}_temp_merged.txt -s
   head -1 ${FILE_PREFIX}_temp_count.txt | sed 's/\t/\n/g' | tail -n+2 > ${FILE_PREFIX}_temp_sample.txt
   fusion_gene_exp_ratio.pl ${snakemake_output[ref_gtf]} ${FILE_PREFIX}_temp_merged.txt ${FILE_PREFIX}_temp_sample.txt ${FILE_PREFIX}
   fusion_gene_monoexonic_finder.pl ${snakemake_output[ref_im_trns_txt]} ${snakemake_output[fusion_mono]}
   echo "gene_id"$'\t'"monoexonic" > ${FILE_PREFIX}_temp_id.txt
   cat ${snakemake_output[fusion_mono]} | sed 's/$/\t1/' >> ${FILE_PREFIX}_temp_id.txt
   merge2tables.pl -t1 ${snakemake_output[fusion_ratio]} -c1 0 -t2 ${FILE_PREFIX}_temp_id.txt -c2 0 -o ${FILE_PREFIX}_temp_merged.txt -s
   sed 's/\t$/\t0/' ${FILE_PREFIX}_temp_merged.txt > ${snakemake_output[fusion_ratio]}
   cp ${FILE_PREFIX}*_fusion_gene_ratio.txt ${FILE_PREFIX}_RECLOCUS_DEBUGOUTPUT_fusion_gene_ratio.txt
   reconstructed_locus.R ${FILE_PREFIX} ${snakemake_input[trackgroups]} 
   cat ${FILE_PREFIX}_reclocus_*_gene_id.txt | sort | uniq > ${snakemake_output[reclocus_id]}
   gtf_reclocus.pl ${snakemake_output[ref_gtf]} ${snakemake_params[refgenome_isop_gtf]} ${snakemake_output[reclocus_id]} ${snakemake_output[reclocus_gtf]}
   gtf_reclocus.pl ${snakemake_output[ref_im_gtf]} ${snakemake_params[refgenome_isop_gtf]} ${snakemake_output[reclocus_id]} ${snakemake_output[reclocus_im_gtf]}

   # Tabulation   
   for attribute in `echo gene_name gene_type status asm_gene_id ref_transcript_id ref_gene_id ref_gene_name ref_gene_type reclocus`; do echo $attribute; done > ${FILE_PREFIX}_temp_attribute_list.txt
   for tag in `echo reconstructed_locus`; do echo $tag; done > ${FILE_PREFIX}_temp_tag_list.txt
   gtf2summary.pl -i ${snakemake_output[reclocus_gtf]} -o ${FILE_PREFIX}_reference_reclocus -a ${FILE_PREFIX}_temp_attribute_list.txt -d ${FILE_PREFIX}_temp_tag_list.txt -g ${snakemake_params[refgenome_fasta]} -j ${snakemake_params[intron_coverage]} -r -t ${snakemake[threads]}
   gtf2summary.pl -i ${snakemake_output[reclocus_im_gtf]} -o ${FILE_PREFIX}_reference_im_reclocus -a ${FILE_PREFIX}_temp_attribute_list.txt -d ${FILE_PREFIX}_temp_tag_list.txt -g ${snakemake_params[refgenome_fasta]} -j ${snakemake_params[intron_coverage]} -r -t ${snakemake[threads]}

   # Select the best prediction from GMST result
   cp ${snakemake_output[reclocus_gtf]} ${FILE_PREFIX}_RECLOCUS_DEBUGOUTPUT_gmst_fnn.gtf
   cp ${snakemake_input[isop_fasta]}     ${FILE_PREFIX}_RECLOCUS_DEBUGOUTPUT_isop_fasta.fasta
   cp ${snakemake_input[gmst_fnn]}     ${FILE_PREFIX}_RECLOCUS_DEBUGOUTPUT_gmst_fnn.fnn
   fasta-reflow.pl  ${snakemake_input[gmst_fnn]} > ${FILE_PREFIX}_gmst.fnn
   GMST2gtf.pl -g ${snakemake_output[reclocus_gtf]} -r ${FILE_PREFIX}_gmst.fnn -f ${snakemake_input[isop_fasta]} -o ${FILE_PREFIX}_reference_reclocus_GMST_CDS -n ${snakemake_params[nmdj_distance]} -e ${snakemake_params[tis_efficiency]}

   # Select the best prediction from CPAT result
   CPAT2gtf.pl -g ${snakemake_output[reclocus_gtf]} -f ${snakemake_input[isop_fasta]} -p ${snakemake_input[cpat_leng]} -s ${snakemake_input[cpat_leng_seqs]} -o ${FILE_PREFIX}_reference_reclocus_CPAT_CDS -n ${snakemake_params[nmdj_distance]} -e ${snakemake_params[tis_efficiency]}

   # Select the best prediction from Transdecoder
   transdecoder2gtf.pl -g ${snakemake_output[reclocus_gtf]} -f ${snakemake_input[isop_fasta]} -r ${snakemake_input[transdecoder]} -o ${FILE_PREFIX}_reference_reclocus_transdecoder_CDS -n ${snakemake_params[nmdj_distance]} -e ${snakemake_params[tis_efficiency]}

   # Select the best prediction among tools
   echo "${FILE_PREFIX}_reference_reclocus_GMST_CDS"$'\t'"GMST" > ${FILE_PREFIX}_temp_prediction_list.txt
   echo "${FILE_PREFIX}_reference_reclocus_CPAT_CDS"$'\t'"CPAT" >> ${FILE_PREFIX}_temp_prediction_list.txt
   echo "${FILE_PREFIX}_reference_reclocus_transdecoder_CDS"$'\t'"Transdecoder" >> ${FILE_PREFIX}_temp_prediction_list.txt
   select_CDS_prediction.pl -i ${FILE_PREFIX}_temp_prediction_list.txt -o ${FILE_PREFIX}_temp_reference_reclocus_CDS.gtf -p ${snakemake_params[min_prob_cpat3]}
   cat ${snakemake_output[reclocus_gtf]} ${FILE_PREFIX}_temp_reference_reclocus_CDS.gtf | sort -k1,1 -k4,4n > ${snakemake_output[reclocus_cds_gtf]}

   # Tabulation   
   for attribute in `echo gene_name gene_type status asm_gene_id ref_gene_id ref_transcript_id ref_gene_name ref_gene_type reclocus cds_source cds_support cds_diversity tis_efficiency`; do echo $attribute; done > ${FILE_PREFIX}_temp_attribute_list.txt
   for tag in `echo reconstructed_locus`; do echo $tag; done > ${FILE_PREFIX}_temp_tag_list.txt
   gtf2summary.pl -i ${snakemake_output[reclocus_cds_gtf]} -o ${FILE_PREFIX}_reference_reclocus_CDS -a ${FILE_PREFIX}_temp_attribute_list.txt -d ${FILE_PREFIX}_temp_tag_list.txt -n ${snakemake_params[nmdj_distance]} -s -g ${snakemake_params[refgenome_fasta]} -j ${snakemake_params[intron_coverage]} -r -t ${snakemake[threads]}

   # Get amino acid sequences
   gtf2fasta_CDS.pl -i ${snakemake_output[reclocus_cds_gtf]} -o ${FILE_PREFIX}_reference_reclocus_CDS -g ${snakemake_params[refgenome_fasta]} -t ${snakemake[threads]}
   
   # Cleanup
   rm -f ${FILE_PREFIX}_temp_attribute_list.txt ${FILE_PREFIX}_temp_tag_list.txt ${FILE_PREFIX}_temp_count.txt ${FILE_PREFIX}_temp_merged.txt ${FILE_PREFIX}_temp_sample.txt ${FILE_PREFIX}_temp_id.txt ${FILE_PREFIX}_temp_prediction_list.txt ${FILE_PREFIX}_temp_reference_reclocus_CDS.gtf

} 2>&1 | tee -a ${snakemake_log[out]}


#----------#
# Run ASEF #
#----------#
{

   # First exon of multiexonic isoforms
   ASEF_FE.pl \
      -q ${snakemake_output[reclocus_cds_gtf]} \
      -r ${snakemake_params[refgenome_isop_gtf]} \
      -d ${snakemake_params[promoter_width]} \
      -o ${FILE_PREFIX}_reference_reclocus_CDS_FE \
      -e ${snakemake_input[tss_bed]}
   
   # Last exon of multiexonic isoforms
   ASEF_LE.pl \
      -q ${snakemake_output[reclocus_cds_gtf]} \
      -r ${snakemake_params[refgenome_isop_gtf]} \
      -o ${FILE_PREFIX}_reference_reclocus_CDS_LE \
      -e ${snakemake_input[tts_bed]}
   
   # Novel exon of multiexonic isoforms
   ASEF_NE.pl \
      -q ${snakemake_output[reclocus_cds_gtf]} \
      -r ${snakemake_params[refgenome_isop_gtf]} \
      -j ${snakemake_params[intron_coverage]} \
      -i \
      -o ${FILE_PREFIX}_reference_reclocus_CDS_NE \
      -t ${snakemake[threads]}
   
   # Novel coding exons
   ASEF_NCE.pl \
      -q ${snakemake_output[reclocus_cds_gtf]} \
      -p \
      -r ${snakemake_params[refgenome_isop_gtf]} \
      -o ${FILE_PREFIX}_reference_reclocus_CDS_NCE \
      -t ${snakemake[threads]}
   
   # Comparison between novel exons and novel coding exons
   ASEF_NCE_NE_comparison.pl \
      -c ${snakemake_output[reclocus_NCE_gtf]} \
      -e ${snakemake_output[reclocus_NE_gtf]} \
      -o ${snakemake_output[reclocus_NECDS]}

} 2>&1 | tee -a ${snakemake_log[out]}



#-------------------------#
# Convert classifications #
#-------------------------#
{
   
   # Merge isoPropeller and ASEF tables
   merge2tables.pl \
      -t1 ${snakemake_output[reclocus_cds_trns_txt]} \
      -c1 0 \
      -t2 ${FILE_PREFIX}_reference_reclocus_CDS_NE_transcript.txt \
      -c2 0 \
      -o ${FILE_PREFIX}_temp_merged.txt \
      -s
   
   # Convert classifications
   classification_conversion.pl \
      ${FILE_PREFIX}_temp_merged.txt \
      ${snakemake_output[reclocus_refined]}

   # Cleanup
   rm -f ${FILE_PREFIX}_temp_merged.txt \

}
