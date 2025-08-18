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
   mkdir -p ${OUT_FOLDER}

} 2>&1 | tee -a ${snakemake_log[out]}


{
   #------------------------------------------#
   # ADD EXTRA ATTRIBUTES TO PATCHED GTF FILE #
   #------------------------------------------#
   
   # Add/replace gene_id
   cut-by-ids -f ${snakemake_input[reclocus_cds_trns_txt]} transcript_id gene_id  | perl -pe 's/\|/:/g' > ${OUT_FOLDER}/map_gene_id.tmp
   gtf-addreplace-attributes.pl -a transcript_id -m ${OUT_FOLDER}/map_gene_id.tmp   -r gene_id            ${snakemake_input[reclocus_cds_gtf]}  > ${FILE_PREFIX}_tmp1.gtf
   rm -f ${OUT_FOLDER}/map_gene_id.tmp
   
   # Add/replace gene_name
   cut-by-ids -f ${snakemake_input[reclocus_cds_trns_txt]} transcript_id gene_name | perl -pe 's/\|/:/g' > ${OUT_FOLDER}/map_gene_name.tmp
   gtf-addreplace-attributes.pl -a transcript_id -m ${OUT_FOLDER}/map_gene_name.tmp -r gene_name                       ${FILE_PREFIX}_tmp1.gtf  > ${FILE_PREFIX}_tmp2.gtf
   rm -f ${OUT_FOLDER}/map_gene_name.tmp ${FILE_PREFIX}_tmp1.gtf
   
   # Add/replace gene_biotype
   cut-by-ids -f ${snakemake_input[reclocus_cds_trns_txt]} transcript_id gene_type | perl -pe 's/\|/:/g' > ${OUT_FOLDER}/map_gene_type.tmp
   gtf-addreplace-attributes.pl -a transcript_id -m ${OUT_FOLDER}/map_gene_type.tmp -r gene_biotype                    ${FILE_PREFIX}_tmp2.gtf  > ${FILE_PREFIX}_tmp3.gtf
   rm -f ${OUT_FOLDER}/map_gene_type.tmp  ${FILE_PREFIX}_tmp2.gtf 

   # Add niap_structural_category
   cut-by-ids -f ${snakemake_input[reclocus_refined]} 'transcript_id' niap_structural_category > ${OUT_FOLDER}/map_niap_cat.tmp
   gtf-addreplace-attributes.pl -a transcript_id -m ${OUT_FOLDER}/map_niap_cat.tmp  -r niap_structural_category        ${FILE_PREFIX}_tmp3.gtf  > ${FILE_PREFIX}_tmp4.gtf
   rm -f ${OUT_FOLDER}/map_niap_cat.tmp  ${FILE_PREFIX}_tmp3.gtf
   
   # Add niap_structural_subcategory
   cut-by-ids -f ${snakemake_input[reclocus_refined]} 'transcript_id' niap_structural_subcategory > ${OUT_FOLDER}/map_niap_subcat.tmp
   gtf-addreplace-attributes.pl -a transcript_id -m ${OUT_FOLDER}/map_niap_subcat.tmp  -r niap_structural_subcategory  ${FILE_PREFIX}_tmp4.gtf  > ${FILE_PREFIX}_tmp5.gtf
   rm -f ${OUT_FOLDER}/map_niap_subcat.tmp ${FILE_PREFIX}_tmp4.gtf
   
   # Add cds_type
   cut-by-ids -f ${snakemake_input[reclocus_cds_trns_txt]} transcript_id cds_type > ${OUT_FOLDER}/map_cds_type.tmp
   gtf-addreplace-attributes.pl -a transcript_id -m ${OUT_FOLDER}/map_cds_type.tmp  -r cds_type                        ${FILE_PREFIX}_tmp5.gtf  > ${FILE_PREFIX}_tmp6.gtf
   rm -f ${OUT_FOLDER}/map_cds_type.tmp ${FILE_PREFIX}_tmp5.gtf
   
   # Add lines for gene entries (requires gffread)
   gtf-add-gene-entries.pl ${FILE_PREFIX}_tmp6.gtf > ${snakemake_output[gtf_final]}
   rm -f ${FILE_PREFIX}_tmp6.gtf
   
   #---------------#
   # ADD EXON RANK #
   #---------------#
   
   # Give everything a rank ID so we can resort in the same order later
   awkt '{print $0, out++}' ${snakemake_output[gtf_final]} > ${FILE_PREFIX}_tmp7.gtf.ranked
   awkt '$3=="exon"' ${FILE_PREFIX}_tmp7.gtf.ranked | perl -pe 's/transcript_id "//; s/";.*\t/\t/' > ${FILE_PREFIX}_tmp7.gtf.ranked.stripped
   
   # Get ranking for pos and neg stranded exons
   awkt '$7=="+"' ${FILE_PREFIX}_tmp7.gtf.ranked.stripped \
      | sortt -k1,1 -k4,4n \
      | awkt '{out[$9]++; print $10,out[$9]}' \
      > ${FILE_PREFIX}_tmp7_exon_rank.txt
   
   awkt '$7=="-"' ${FILE_PREFIX}_tmp7.gtf.ranked.stripped \
      | sortt -k1,1 -k4,4nr \
      | awkt '{out[$9]++; print $10,out[$9]}' \
      >> ${FILE_PREFIX}_tmp7_exon_rank.txt
   
   # Put it all back together again
   join-by-ids -a 1 -1 10 -2 1 ${FILE_PREFIX}_tmp7.gtf.ranked ${FILE_PREFIX}_tmp7_exon_rank.txt \
      | sortt -k 10,10n \
      | awkt '{print $1,$2,$3,$4,$5,$6,$7,$8,$9 " exon_number " $11 ";"}' \
      > ${FILE_PREFIX}_tmp8.gtf
   perl -pi -e 's/ exon_number ;$//' ${FILE_PREFIX}_tmp8.gtf
   
   
   #-----------------#
   # FIX STOP CODONS #
   #-----------------#
   batch-gtf-remove-stopcodon-from-CDS.sh ${FILE_PREFIX}_tmp8.gtf ${FILE_PREFIX}_patched_extra_stopfix.gtf
   rm -f ${FILE_PREFIX}_tmp7* ${FILE_PREFIX}_tmp8.gtf
   
   
   #---------------------------------#
   # GENERATE PER-SAMPLE TRACK FILES #
   #---------------------------------#
   for track in `head -n1 ${snakemake_input[iso_count_matrix]} | cut --complement -f 1`
   do
      trackbase="${FILE_PREFIX}_${track/flcount_/sample_}"
      cut-by-ids -f ${snakemake_input[iso_count_matrix]} TranscriptID $track | awkt '$2>0 && $1 !~ /TranscriptID/ {print $1, $2}' > ${trackbase}.ids
      gtf-filter-attributes.pl -a transcript_id -m ${trackbase}.ids ${snakemake_output[gtf_final]} > ${trackbase}.gtf
      gtf2gff.pl -a transcript_id,gene_name,gene_biotype,cds_type -f exon,CDS ${trackbase}.gtf > ${trackbase}.gff
      gff2bed.pl ${trackbase}.gff > ${trackbase}.bed
      awkt '{split($4,a,"|"); print a[1], $0}' ${trackbase}.bed > ${trackbase}.bed.score1
      join-by-ids ${trackbase}.ids ${trackbase}.bed.score1 | awkt '{$7=$2; print $0}' | cut --complement -f1,2 | sortt -k1,1 -k2,2n -k3,3n > ${trackbase}.bed
      rm -f ${trackbase}.ids ${trackbase}.gff ${trackbase}.bed.score1
      gzip ${trackbase}.gtf
   done

   #--------------------------------#
   # GENERATE PER-GROUP TRACK FILES #
   #--------------------------------#
   for group in `cut -f2 ${snakemake_input[trackgroups]} | sort | uniq | xargs`
   do
      trackbase="${FILE_PREFIX}_${group}"
      awkt -v group=${group} '$2==group{print $1}' ${snakemake_input[trackgroups]} > ${trackbase}.countfilt
      intersect-cols-by-ids -ff ${snakemake_input[iso_count_matrix]} -if ${trackbase}.countfilt -fo 2 \
         | awkt '{sum=0; for(i=2;i<=NF;i++){sum+=$i}; print $1, sum }' \
         | awkt '$2>0 && $1 !~ /TranscriptID/ {print $1, $2}' > ${trackbase}.ids
      gtf-filter-attributes.pl -a transcript_id -m ${trackbase}.ids ${snakemake_output[gtf_final]} > ${trackbase}.gtf
      gtf2gff.pl -a transcript_id,gene_name,gene_biotype,cds_type -f exon,CDS ${trackbase}.gtf > ${trackbase}.gff
      gff2bed.pl ${trackbase}.gff > ${trackbase}.bed
      awkt '{split($4,a,"|"); print a[1], $0}' ${trackbase}.bed > ${trackbase}.bed.score1
      join-by-ids ${trackbase}.ids ${trackbase}.bed.score1 | awkt '{$7=$2; print $0}' | cut --complement -f1,2 | sortt -k1,1 -k2,2n -k3,3n > ${trackbase}.bed
      rm -f ${trackbase}.ids ${trackbase}.gff ${trackbase}.bed.score1 ${trackbase}.countfilt
      gzip ${trackbase}.gtf
   done
   
} 2>&1 | tee -a ${snakemake_log[out]}

