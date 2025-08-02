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
OUT_FOLDER="${snakemake_params[output_dir]}"
FILE_PREFIX="${OUT_FOLDER}/${snakemake_params[prefix]}"


########
# MAIN #
########

#----------------#
# Prepare inputs #
#----------------#
{
   echo -e "\n ## Prepare inputs ## \n"
   
   # Set up output folder
   rm -rf ${OUT_FOLDER}
   mkdir ${OUT_FOLDER}
   
   # Convert sqanti3 gtf file to bed format
   gtf2gff.pl -a transcript_id ${snakemake_input[sq3_gtf]}  > ${FILE_PREFIX}_corrected.tmp.gff
   gff2bed.pl ${FILE_PREFIX}_corrected.tmp.gff              > ${FILE_PREFIX}_corrected.tmp.bed
   
   # Collapse all reference genome transcripts based on any exon overlap on the same strand
   gtf2gff.pl -a gene_id,transcript_id ${snakemake_params[refgenome_gtf]} > ${FILE_PREFIX}_reference-annotations.gff   
   gff2bed.pl ${FILE_PREFIX}_reference-annotations.gff                    > ${FILE_PREFIX}_reference-annotations.bed
   collapseBed.pl ${FILE_PREFIX}_reference-annotations.bed                > ${FILE_PREFIX}_reference-annotations_collapsed.bed

} 2>&1 | tee -a ${snakemake_log[out]}

#------------------------------#
# Get corrected cage distances #
#------------------------------#
{
   echo -e "\n ## Get corrected cage distances ## \n"
   
   # Get additional input files for bedtools
   awkt '$6 == "+" {$3=$2; print $0} $6 == "-" {$2=$3; print $0}' ${FILE_PREFIX}_corrected.tmp.bed | sortt -k1,1 -k2,2n | uniq | cut -f 1-6 | grep -P '^chr' > ${FILE_PREFIX}_closest_ISOF_5prime.bed
   awkt '$6 == "+" {$2=$3; print $0} $6 == "-" {$3=$2; print $0}' ${FILE_PREFIX}_corrected.tmp.bed | sortt -k1,1 -k2,2n | uniq | cut -f 1-6 | grep -P '^chr' > ${FILE_PREFIX}_closest_ISOF_3prime.bed
   sortt -k1,1 -k2,2n ${snakemake_params[cage_peaks_allmerge]} | awkt '{$4=out++; print $0}' | cut -f 1-6 > ${FILE_PREFIX}_closest_CAGE.bed
   sortt -k1,1 -k2,2n ${snakemake_params[polya_bed]}           | awkt '{$4=out++; print $0}' | cut -f 1-6 > ${FILE_PREFIX}_closest_POLY.bed
   
   # Get distances from 5' to closest cage peak and 3' to closest polyA peak
   bedtools closest -a ${FILE_PREFIX}_closest_ISOF_5prime.bed -b ${FILE_PREFIX}_closest_CAGE.bed -t first -s -D a -nonamecheck | cut -f 4,13 > ${OUT_FOLDER}/bedtools_closest_cage_peak
   add-header ${OUT_FOLDER}/bedtools_closest_cage_peak "#PBID" closest_cage_peak
   bedtools closest -a ${FILE_PREFIX}_closest_ISOF_3prime.bed -b ${FILE_PREFIX}_closest_POLY.bed -t first -s -D a -nonamecheck | cut -f 4,13 > ${OUT_FOLDER}/bedtools_closest_polyA_peak
   add-header ${OUT_FOLDER}/bedtools_closest_polyA_peak "#PBID" closest_polyA_peak
   
   # Join it all together
   multi-join-by-ids -e "NA" ${OUT_FOLDER}/bedtools_closest_cage_peak ${OUT_FOLDER}/bedtools_closest_polyA_peak \
      | perl -pe 's/^#joinID/isoform/' \
      > ${FILE_PREFIX}.renamed_corrected.closestPeaks.txt
   
   # Cleanup intermediary files for this module
   rm -f ${FILE_PREFIX}_bedtools_closest_cage_peak \
         ${FILE_PREFIX}_bedtools_closest_polyA_peak \
         ${FILE_PREFIX}_closest_ISOF_5prime.bed \
         ${FILE_PREFIX}_closest_ISOF_3prime.bed \
         ${FILE_PREFIX}_closest_CAGE.bed \
         ${FILE_PREFIX}_closest_POLY.bed \
         ${OUT_FOLDER}/bedtools_closest_cage_peak \
         ${OUT_FOLDER}/bedtools_closest_polyA_peak
         
} 2>&1 | tee -a ${snakemake_log[out]}


##---------------------------#
## Get extension information #
##---------------------------#
{
   echo -e "\n ## Get extension information ## \n"
   
   # Run the extension script to find any isoforms that extend a known gene
   find-gene-joins-and-extensions.pl \
      -i ${FILE_PREFIX}_corrected.tmp.bed \
      -r ${FILE_PREFIX}_reference-annotations_collapsed.bed \
      -a ${snakemake_params[refgenome_anno]} \
      > ${FILE_PREFIX}.renamed_corrected.extensions.txt
   add-header ${FILE_PREFIX}.renamed_corrected.extensions.txt isoform extension_type extension_gene_id extension_gene_name extension_gene_typeA extension_gene_typeB

} 2>&1 | tee -a ${snakemake_log[out]}


##------------------------------#
## Get collapsed novelGene loci #
##------------------------------#
{
   echo -e "\n ## Get collapsed novelGene loci ## \n"

   # Prepare a bed file with all novelGene isoforms
   awkt '$7 ~ /^novelGene/ {print $1}' ${snakemake_input[sq3_classification]} > ${FILE_PREFIX}.collapse-tmp.PBIDs
   intersect-by-ids -ff ${FILE_PREFIX}_corrected.tmp.bed -fc 4 -if ${FILE_PREFIX}.collapse-tmp.PBIDs > ${FILE_PREFIX}.collapse-tmp.PBIDs.bed

   # Collapse the novelgenes by any exon overlap on the same strand, then renumber the outputs to a new set of novelGene IDs
   collapseBed.pl ${FILE_PREFIX}.collapse-tmp.PBIDs.bed \
      | cut -f 4 \
      | awkt '{out++; split($1,a,"|"); for(var in a){print a[var], "novelGene" out}}' \
      > ${FILE_PREFIX}.renamed_corrected.novelGene.txt
   
   # Cleanup intermediary files for this module
   rm -f ${FILE_PREFIX}.collapse-tmp.PBIDs \
         ${FILE_PREFIX}.collapse-tmp.PBIDs.bed

} 2>&1 | tee -a ${snakemake_log[out]}


## -------------------------------------------------#
## Get gene mappings based on best-match transcript #
##--------------------------------------------------#
{
   echo -e "\n ## Get gene mappings based on best-match transcript ## \n"

   # Here we get transcript ID to gene ID mappings from the reference genome
   gtf-count-attributes.pl -a transcript_id,gene_id ${snakemake_params[refgenome_gtf]} \
      | cut -f 1-2 \
      > ${FILE_PREFIX}.renamed_corrected.bestHitFromTranscript.txt

} 2>&1 | tee -a ${snakemake_log[out]}


##-----------------------------------------------------------------------------------------------------------------#
## Update sqanti file structural classifications, extensions, closest peaks, best hits, and add novelgene overlaps #
##-----------------------------------------------------------------------------------------------------------------#
{
   echo -e "\n ## Update sqanti file structural classifications, extensions, closest peaks, best hits, and add novelgene overlaps ## \n"

   # Integrate all information we just collected in a new corrected sqanti output file
   isoseq_sqanti-reclassify.pl \
      -f ${snakemake_input[sq3_classification]} \
      -m ${FILE_PREFIX}.renamed_corrected.novelGene.txt \
      -e ${FILE_PREFIX}.renamed_corrected.extensions.txt \
      -c ${FILE_PREFIX}.renamed_corrected.closestPeaks.txt \
      -t ${FILE_PREFIX}.renamed_corrected.bestHitFromTranscript.txt \
      2> ${snakemake_output[stats_reclassify]} \
      | awkt '$2~/^chr/ && $2!~/^chrM/' \
      | perl -pe 's/^isoform\t/#isoform\t/' \
      > "${snakemake_output[sq3_classpatched]}.tmp1"
   
   # Reclassify novelgene overlaps
   isoseq_sqanti-reclassify_novelgene_overlaps.pl \
      -b \
      -i ${FILE_PREFIX}_corrected.tmp.bed \
      -a "${snakemake_output[sq3_classpatched]}.tmp1" \
      2> ${snakemake_output[stats_novelgene]} \
      > "${snakemake_output[sq3_classpatched]}.tmp2"
      
   # Cleanup intermediary files for this module
   rm -f "${snakemake_output[sq3_classpatched]}.tmp1"

} 2>&1 | tee -a ${snakemake_log[out]}


##-----------------------#
## Annotate fusion genes #
##-----------------------#
{
   echo -e "\n ## Annotate fusion genes ## \n"

   awkt '$6!~/^fusion/ && $51!="joins" && $7 ~/_ENSG/ {out="Not called by any and assigned to more than one gene"}  
         $6!~/^fusion/ && $51!="joins" && $7!~/_ENSG/ {out=""} 
         $6 ~/^fusion/ && $51=="joins"                {out="Called by both"} 
         $6 ~/^fusion/ && $51!="joins"                {out="Called by sqanti"} 
         $6!~/^fusion/ && $51=="joins"                {out="Called by extensions"} 
         $1=="#isoform"                               {out="fusion_combined_evidence"} 
         {print $0, out}' \
         "${snakemake_output[sq3_classpatched]}.tmp2" \
         > ${snakemake_output[sq3_classpatched]}

   # Cleanup intermediary files for this module
   rm -f "${snakemake_output[sq3_classpatched]}.tmp2"

} 2>&1 | tee -a ${snakemake_log[out]}


##-----------------------------------------------------------------#
## Generate a new GTF file that incorporates the new gene mappings #
##-----------------------------------------------------------------#
{
   echo -e "\n ## Generate a new GTF file that incorporates the new gene mappings ## \n"
   
   # Get a temporary file with remapped gene IDs for each isoform
   cut-by-ids -f ${snakemake_output[sq3_classpatched]} isoform associated_gene > ${FILE_PREFIX}_gene-id-remapping.tmp

   # Remap the gene IDs and create new corrected gtf file
   gtf-addreplace-attributes.pl \
      -a transcript_id \
      -r gene_id \
      -m ${FILE_PREFIX}_gene-id-remapping.tmp \
      ${snakemake_input[sq3_gtf]} \
      >  ${FILE_PREFIX}_gene-id-remapping.tmp.gtf
      
   # NIAP and ASEF need proper frame information in the GTF file, so we pass it through gffread
   gffread -T ${FILE_PREFIX}_gene-id-remapping.tmp.gtf > ${snakemake_output[gtf_classpatched]}

   # Summarize the patched GTF file
   gtf2summary.pl -i ${snakemake_output[gtf_classpatched]} -o ${FILE_PREFIX}.txt -s -n 50 -g ${snakemake_params[refgenome_fasta]} -j ${snakemake_params[intron_coverage]} -r -t ${snakemake[threads]}

   # Also create outputs in bed and gff formats
   gtf2gff.pl ${snakemake_output[gtf_classpatched]} > ${snakemake_output[gff_classpatched]}
   gff2bed.pl ${snakemake_output[gff_classpatched]} > ${snakemake_output[bed_classpatched]}

   # Cleanup intermediary files for this module
   rm -f ${FILE_PREFIX}_gene-id-remapping.tmp \
         ${FILE_PREFIX}_gene-id-remapping.tmp.gtf

} 2>&1 | tee -a ${snakemake_log[out]}


##------------------------------------------------------------------#
## Final cleanup of intermediary files used by more than one module #
##------------------------------------------------------------------#
{
   echo -e "\n ## Final cleanup of intermediary files ## \n"

   rm -f ${FILE_PREFIX}_corrected.tmp.gff \
         ${FILE_PREFIX}_corrected.tmp.bed \
         ${FILE_PREFIX}_reference-annotations.gff \
         ${FILE_PREFIX}_reference-annotations.bed \
         ${FILE_PREFIX}_reference-annotations_collapsed.bed \
         ${FILE_PREFIX}.renamed_corrected.bestHitFromTranscript.txt \
         ${FILE_PREFIX}.renamed_corrected.closestPeaks.txt \
         ${FILE_PREFIX}.renamed_corrected.extensions.txt \
         ${FILE_PREFIX}.renamed_corrected.novelGene.txt

} 2>&1 | tee -a ${snakemake_log[out]}
