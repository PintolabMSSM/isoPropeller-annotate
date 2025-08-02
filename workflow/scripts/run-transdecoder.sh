#!/usr/bin/env bash

###############
# ENVIRONMENT #
###############

# Set paths from snakemake parameters
export PATH="${snakemake_params[path]}:$PATH"
export PERL5LIB="${snakemake_params[perl5lib]}:$PERL5LIB"

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
   rm -rf ${OUT_FOLDER}
   mkdir -p ${OUT_FOLDER}

} 2>&1 | tee -a ${snakemake_log[out]}


#-----------------------------------#
# Run transdecoder with pfam search #
#-----------------------------------#
{
   # Change into output dir
   cd ${OUT_FOLDER}
   
   # Find longest ORFs on the sense strand
   TransDecoder.LongOrfs -m ${snakemake_params[minorf_transdecoder]} -S -t ../${snakemake_input[sq3_fasta]} --output_dir ./
   
   # Split longest ORF peptide file into 24 chunks for parallelization
   SPLITLIST="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24"
   fasta-splitter-random.sh "${snakemake_params[prefix]}_corrected.fasta.transdecoder_dir/longest_orfs.pep" 24 split
   
   # Run a parallelized hmmscan search with pfam domains
   parallel --jobs 24 hmmsearch -E 1e-10 --domtblout split_{}.domtblout ${snakemake_params[pfam_a_hmm]} split_{}.fa ::: ${SPLITLIST}
   cat split_*.domtblout > "${snakemake_params[prefix]}_corrected.fasta.transdecoder.domtblout"
   
   # Run a parallelized blastp search against uniprot_sprot using blast+. This option is more sensitive but takes about 100x longer
   # parallel --jobs 24  blastp -query "${snakemake_params[prefix]}_corrected.fasta.transdecoder_dir/longest_orfs.pep" -db ${snakemake_params[uniprot_sprot]} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 2 -out split_{}.blastp ::: ${SPLITLIST}   
   # cat split_*.domtblout > "${snakemake_params[prefix]}_corrected.fasta.transdecoder.blastp"
   
   # Run parallelized blastp search against uniprot_sport using diamondblast. Much faster and allows use of uniref90 database
   diamond blastp --sensitive --query "${snakemake_params[prefix]}_corrected.fasta.transdecoder_dir/longest_orfs.pep" --db ${snakemake_params[uniref90]} --max-target-seqs 1 --outfmt 6 --evalue 1e-5 --threads ${snakemake[threads]} --out "${snakemake_params[prefix]}_corrected.fasta.transdecoder.blastp"
   
   # Run predict and retain pfam and/or blast hits
   TransDecoder.Predict --retain_pfam_hits "${snakemake_params[prefix]}_corrected.fasta.transdecoder.domtblout" --retain_blastp_hits "${snakemake_params[prefix]}_corrected.fasta.transdecoder.blastp" -t ../${snakemake_input[sq3_fasta]} -O ./
   
   # Let Snakemake know that we're done here
   touch "${snakemake_params[prefix]}_transdecoder.check"

   # Cleanup
   rm -f split_*.fa split_*.domtblout
   
} 2>&1 | tee -a ${snakemake_log[out]}

