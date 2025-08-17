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

#----------------------------------------------#
# Prepare output folder and files for overlaps #
#----------------------------------------------#
{
   # Set up output folder
   mkdir ${OUT_FOLDER}
   
   # Prep data for exon-level overlaps
   gtf2gff.pl -a transcript_id ${snakemake_input[reclocus_cds_gtf]}         > ${FILE_PREFIX}_corrected.tmp-exon.gff
   gff2bed.pl ${FILE_PREFIX}_corrected.tmp-exon.gff                         > ${FILE_PREFIX}_corrected.tmp-exon.bed
   
   # Prep data for CDS-level overlaps
   gtf2gff.pl -f CDS -a transcript_id ${snakemake_input[reclocus_cds_gtf]}  > ${FILE_PREFIX}_corrected.tmp-CDS.gff
   gff2bed.pl ${FILE_PREFIX}_corrected.tmp-CDS.gff                          > ${FILE_PREFIX}_corrected.tmp-CDS.bed

} 2>&1 | tee -a ${snakemake_log[out]}


#------------------------------------#
# Get overlaps with genomic elements #
#------------------------------------#
{

   # Segmental duplications
   bedtools intersect -nonamecheck -wao -split -a ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${snakemake_params[segdups_bed]} \
      | awkt '{out[$4]+=$NF} END{print "#pbid\toverlap_segdups_exon_bp"; for(var in out){print var, out[var]}}' \
      > ${OUT_FOLDER}/overlap_segdups_exon
   
   bedtools intersect -nonamecheck -wao -split -a ${FILE_PREFIX}_corrected.tmp-CDS.bed -b ${snakemake_params[segdups_bed]} \
      | awkt '{out[$4]+=$NF} END{print "#pbid\toverlap_segdups_CDS_bp"; for(var in out){print var, out[var]}}' \
      > ${OUT_FOLDER}/overlap_segdups_CDS
   
   # Repeatmasker
   bedtools intersect -nonamecheck -split -a ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${snakemake_params[rmsk_bed]} -wao \
      | awkt '{out[$4] += $19} END{print "#pbid\toverlap_repeatmasker_exon_bp"; for(var in out) {print var, out[var]}}'  \
      > ${OUT_FOLDER}/overlap_repeatmasker_exon
   
   bedtools intersect -nonamecheck -split -a ${FILE_PREFIX}_corrected.tmp-CDS.bed -b ${snakemake_params[rmsk_bed]} -wao \
      | awkt '{out[$4] += $19} END{print "#pbid\toverlap_repeatmasker_CDS_bp"; for(var in out) {print var, out[var]}}'  \
      > ${OUT_FOLDER}/overlap_repeatmasker_CDS
   
   # SINE, LINE and LTR element overlaps
   bedtools intersect -nonamecheck -split -a ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${snakemake_params[rmsk_sel_bed]} -wao \
      | awkt '   {out[$4]++; sine[$4] += $16~/SINE/ ? $19 : 0; line[$4] += $16~/LINE/ ? $19 : 0; ltr[$4]  += $16~/LTR/  ? $19 : 0; dna[$4]  += $16~/DNA/  ? $19 : 0;} END{print "#pbid\toverlap_SN|LN|LTR|DNA_exon_bp"; for(i in out){print i, sine[i]"|"line[i]"|"ltr[i]"|"dna[i]}}' \
      > ${OUT_FOLDER}/overlap_SINE-LINE-LTR_exon
   
   bedtools intersect -nonamecheck -split -a ${FILE_PREFIX}_corrected.tmp-CDS.bed -b ${snakemake_params[rmsk_sel_bed]} -wao \
      | awkt '   {out[$4]++; sine[$4] += $16~/SINE/ ? $19 : 0; line[$4] += $16~/LINE/ ? $19 : 0; ltr[$4]  += $16~/LTR/  ? $19 : 0; dna[$4]  += $16~/DNA/  ? $19 : 0;} END{print "#pbid\toverlap_SN|LN|LTR|DNA_CDS_bp"; for(i in out){print i, sine[i]"|"line[i]"|"ltr[i]"|"dna[i]}}' \
      > ${OUT_FOLDER}/overlap_SINE-LINE-LTR_CDS
   
   # ultraconserved element overlaps
   bedtools intersect -nonamecheck -s -abam ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${snakemake_params[ultracons]} -wao -split -bed \
      | awkt '{b[$4]+=$19} END{print "#pbid\toverlap_ultraconserved_exon_bp"; for(i in b){print i,b[i]}}' \
      > ${OUT_FOLDER}/overlap_ultraconserved_exon
   
   bedtools intersect -nonamecheck -s -abam ${FILE_PREFIX}_corrected.tmp-CDS.bed -b ${snakemake_params[ultracons]} -wao -split -bed \
      | awkt '{b[$4]+=$19} END{print "#pbid\toverlap_ultraconserved_CDS_bp"; for(i in b){print i,b[i]}}' \
      > ${OUT_FOLDER}/overlap_ultraconserved_CDS
   
   # PHYLOCSF_v31
   bedtools intersect -nonamecheck -s -abam ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${snakemake_params[phylocsf_v31]} -wao -split -bed \
      | awkt '{b[$4]+=$26} END{print "#pbid\toverlap_phyloCSF_novel_v31_exon"; for(i in b){print i,b[i]}}' \
      > ${OUT_FOLDER}/overlap_phyloCSF_novel_v31_exon
   
   bedtools intersect -nonamecheck -s -abam ${FILE_PREFIX}_corrected.tmp-CDS.bed -b ${snakemake_params[phylocsf_v31]} -wao -split -bed \
      | awkt '{b[$4]+=$26} END{print "#pbid\toverlap_phyloCSF_novel_v31_CDS"; for(i in b){print i,b[i]}}' \
      > ${OUT_FOLDER}/overlap_phyloCSF_novel_v31_CDS
   
   # PHYLOCSF_v35
   bedtools intersect -nonamecheck -s -abam ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${snakemake_params[phylocsf_v35]} -wao -split -bed \
      | awkt '{b[$4]+=$26} END{print "#pbid\toverlap_phyloCSF_novel_v35_exon"; for(i in b){print i,b[i]}}' \
      > ${OUT_FOLDER}/overlap_phyloCSF_novel_v35_exon
      
   bedtools intersect -nonamecheck -s -abam ${FILE_PREFIX}_corrected.tmp-CDS.bed -b ${snakemake_params[phylocsf_v35]} -wao -split -bed \
      | awkt '{b[$4]+=$26} END{print "#pbid\toverlap_phyloCSF_novel_v35_CDS"; for(i in b){print i,b[i]}}' \
      > ${OUT_FOLDER}/overlap_phyloCSF_novel_v35_CDS 
   
   # Put it all together
   multi-join-by-ids -e 0 \
      ${OUT_FOLDER}/overlap_segdups_exon \
      ${OUT_FOLDER}/overlap_segdups_CDS \
      ${OUT_FOLDER}/overlap_repeatmasker_exon \
      ${OUT_FOLDER}/overlap_repeatmasker_CDS  \
      ${OUT_FOLDER}/overlap_SINE-LINE-LTR_exon \
      ${OUT_FOLDER}/overlap_SINE-LINE-LTR_CDS \
      ${OUT_FOLDER}/overlap_ultraconserved_exon \
      ${OUT_FOLDER}/overlap_ultraconserved_CDS \
      ${OUT_FOLDER}/overlap_phyloCSF_novel_v31_exon \
      ${OUT_FOLDER}/overlap_phyloCSF_novel_v31_CDS \
      ${OUT_FOLDER}/overlap_phyloCSF_novel_v35_exon \
      ${OUT_FOLDER}/overlap_phyloCSF_novel_v35_CDS \
      > ${FILE_PREFIX}_genomic_element_overlaps.txt
   
   # Cleanup intermediary files
   rm -f ${OUT_FOLDER}/overlap_*_CDS \
         ${OUT_FOLDER}/overlap_*_exon

} 2>&1 | tee -a ${snakemake_log[out]}


#-----------------------------------------------#
# Get overlaps with common structural variation #
#-----------------------------------------------#
{

   # Split datasets by DEL, DUP, OTHER and make non-redundant
   # AF in col 5, EUR_AF in col 6
   for freq in 0.0001 0.005 0.01
   do
      grep "_DUP_"            ${snakemake_params[sv_ctrl]} | awkt -v freq=$freq '$5<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_all_ctrl_DUP_${freq}.bed
      grep "_DEL_"            ${snakemake_params[sv_ctrl]} | awkt -v freq=$freq '$5<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_all_ctrl_DEL_${freq}.bed
      grep -vP "_DEL_|_DUP_"  ${snakemake_params[sv_ctrl]} | awkt -v freq=$freq '$5<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_all_ctrl_OTH_${freq}.bed
      grep "_DUP_"            ${snakemake_params[sv_nneu]} | awkt -v freq=$freq '$5<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_all_nneu_DUP_${freq}.bed
      grep "_DEL_"            ${snakemake_params[sv_nneu]} | awkt -v freq=$freq '$5<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_all_nneu_DEL_${freq}.bed
      grep -vP "_DEL_|_DUP_"  ${snakemake_params[sv_nneu]} | awkt -v freq=$freq '$5<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_all_nneu_OTH_${freq}.bed
      
      grep "_DUP_"            ${snakemake_params[sv_ctrl]} | awkt -v freq=$freq '$6<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_eur_ctrl_DUP_${freq}.bed
      grep "_DEL_"            ${snakemake_params[sv_ctrl]} | awkt -v freq=$freq '$6<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_eur_ctrl_DEL_${freq}.bed
      grep -vP "_DEL_|_DUP_"  ${snakemake_params[sv_ctrl]} | awkt -v freq=$freq '$6<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_eur_ctrl_OTH_${freq}.bed
      grep "_DUP_"            ${snakemake_params[sv_nneu]} | awkt -v freq=$freq '$6<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_eur_nneu_DUP_${freq}.bed
      grep "_DEL_"            ${snakemake_params[sv_nneu]} | awkt -v freq=$freq '$6<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_eur_nneu_DEL_${freq}.bed
      grep -vP "_DEL_|_DUP_"  ${snakemake_params[sv_nneu]} | awkt -v freq=$freq '$6<freq' | cut -f 1-4 | sortt -k1,1 -k2,2n -k3,3n | bedtools merge | awkt '{print $1,$2,$3,"out"count++}' > ${OUT_FOLDER}/SV_eur_nneu_OTH_${freq}.bed 
   done
   
   # Parse all the overlaps
   for i in ${OUT_FOLDER}/SV_???_????_???_*.bed
   do
      name=`echo $i | perl -pe 's/\.bed$//'`
      bedtools merge -c 4 -o collapse -i ${i} > ${name}_nr.bed
      bedtools intersect -a ${FILE_PREFIX}_corrected.tmp-exon.bed -b ${i} -wao -split -bed \
         | awkt -v name=$name '{b[$4]+=$17} END{print "#pbid\t" name; for(i in b){print i,b[i]}}' \
         > ${name}_ovlpbp
      rm -f ${i} ${name}_nr.bed
   done
   
   # Merge data
   multi-join-by-ids -e 0 ${OUT_FOLDER}/SV_???_????_???_*_ovlpbp > ${FILE_PREFIX}_SV_overlaps.txt
   
   # Cleanup intermediary files
   rm -f ${OUT_FOLDER}/SV_???_????_???_*_ovlpbp

} 2>&1 | tee -a ${snakemake_log[out]}


#---------------#
# Final cleanup #
#---------------#
{
   # We can get rid of the CDS/exon files
   rm -f ${FILE_PREFIX}_corrected.tmp-exon.gff \
         ${FILE_PREFIX}_corrected.tmp-exon.bed \
         ${FILE_PREFIX}_corrected.tmp-CDS.gff \
         ${FILE_PREFIX}_corrected.tmp-CDS.bed \

} 2>&1 | tee -a ${snakemake_log[out]}
