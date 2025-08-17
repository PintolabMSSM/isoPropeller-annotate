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

#--------------------------------------------------------------------------------#
# Generate the set of inputs that are needed for isoformswitchanalyzer analysis  #
#--------------------------------------------------------------------------------#
{
   # Set up output folder
   mkdir ${OUT_FOLDER}
   
   # Generate the counts file, the filtered counts file and the matching design and TPM files
   awkt 'NR==1 && /^#/ {for (i=2; i<=NF; i++) if (i==2) printf $i; else printf "\t"$i; printf "\n"} NR>1' ${snakemake_input[iso_count_matrix]} > ${FILE_PREFIX}_exp_counts.tmp
   R --vanilla <<EOF
      # Read count and group/design file
      d     = read.delim('${FILE_PREFIX}_exp_counts.tmp', check.names=F)
      g     = read.delim('${snakemake_input[trackgroups]}', header=F)
      colnames(g) = c("sampleID","condition")

      # Order the group/design and count file by condition, then sampleID
      g     = g[ order(g[,2], g[,1]) ,]
      d     = d[ , g[,1] ]
      
      # Create filtered versions of the count and design files
      #  - Keep only samples with a minimal total read count of 'isoswitch_min_count' (default: 5000)
      #  - Then, keep only groups with at least 2 samples after filtering on 'isoswitch_min_count'
      flt   = colSums(d) > ${snakemake_params[isoswitch_min_count]}
      grpc  = table(g[flt,2])
      grpf  = names(grpc)[grpc>1]
      gf    = g[ flt & g[,2] %in% grpf, ]
      df    = d[ , gf[,1] ]

      # Make TPM abundance matrices
      dtpm  = t(t(as.matrix(d))/colSums(d))*1000000
      dftpm = t(t(as.matrix(df))/colSums(df))*1000000
      
      # Write outputs
      write.table(g,     file="${FILE_PREFIX}_exp_design.txt",          row.names=F, col.names=T, sep="\t", quote=F)
      write.table(gf,    file="${FILE_PREFIX}_exp_design_filt.txt",     row.names=F, col.names=T, sep="\t", quote=F)
      write.table(d,     file="${FILE_PREFIX}_exp_counts.txt",      row.names=T, col.names=T, sep="\t", quote=F)
      write.table(df,    file="${FILE_PREFIX}_exp_counts_filt.txt", row.names=T, col.names=T, sep="\t", quote=F)
      write.table(dtpm,  file="${FILE_PREFIX}_exp_TPM.txt",         row.names=T, col.names=T, sep="\t", quote=F)
      write.table(dftpm, file="${FILE_PREFIX}_exp_TPM_filt.txt",    row.names=T, col.names=T, sep="\t", quote=F)
EOF

   # Gather the stopfix file from the tracks folder : annots.gtf
   cp ${snakemake_input[gtf_stopfix]} ${FILE_PREFIX}_exp_annots.gtf
     
   # Copy cpat3 output and include a reformatted version to cpatv2: cpat2.txt. cpat3.txt
   cp ${snakemake_input[cpat_prob_best]} ${FILE_PREFIX}_exp_cpat3.txt
   echo -e "mRNA_size\tORF_size\tFickett_score\tHexamer_score\tcoding_prob" > ${FILE_PREFIX}_exp_cpat2.txt
   cut-by-ids -f ${FILE_PREFIX}_exp_cpat3.txt seq_ID mRNA ORF Fickett Hexamer Coding_prob \
      |  grep -vP '^seq_ID\t' \
      >> ${FILE_PREFIX}_exp_cpat2.txt
   
   # Copy pfamscan output (if it exists) : isoforms.pfamscan
   if [[ -e ${snakemake_input[pfamscan_output]} ]]
   then
      cp ${snakemake_params[pfamscan_output]} ${FILE_PREFIX}_exp_pfam.pfamscan
   fi
   
   # Copy the isoform faa file : isoforms.faa
   cp ${snakemake_input[reclocus_cds_aa]} ${FILE_PREFIX}_exp_isoforms.faa
   
   # Derive the isoform nucleotide sequence fasta file from the stopfix gtf file : isoforms-nt.fasta
   gtf2gff.pl ${FILE_PREFIX}_exp_annots.gtf > ${FILE_PREFIX}_exp_annots.gff
   gff2bed.pl ${FILE_PREFIX}_exp_annots.gff > ${FILE_PREFIX}_exp_annots.bed
   bedtools getfasta \
      -fi       ${snakemake_params[refgenome_fasta]} \
      -bed      ${FILE_PREFIX}_exp_annots.bed \
      -nameOnly \
      -s        \
      -split    \
      | perl -pe 's/\([+-]\)$//' \
      > ${FILE_PREFIX}_exp_isoforms-nt.fasta

} 2>&1 | tee -a ${snakemake_log[out]}

