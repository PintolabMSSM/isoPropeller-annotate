#-----------------------------------------------------#
# GET INPUT FILES FOR ISOFORMSWITCHANALYZER ANALYSIS  #
#-----------------------------------------------------#
rule get_isoformswitchanalyzer_input_files:
   message: "Gather the set of files that are needed as inputs for isoformswitchanalyzer"
   input:
      iso_count_matrix    = "03_isoform_counts/{prefix}_isoform-counts.txt",
      gtf_stopfix         = "12_tracks/{prefix}_patched_extra_stopfix.gtf",
      trackgroups         = "{prefix}.trackgroups",
      cpat_prob_best      = "05_cpatv3/{prefix}_corrected.cpatv3p18.ORF_prob.best.tsv",
      pfam_launched       = "07_pfam/{prefix}_jobs-submitted.check",
      reclocus_cds_aa     = "09_niap_asef/{prefix}_reference_reclocus_CDS_aa.fa"      
   output:
      exp_counts          = "16_isoformswitchanalyzer_inputs/{prefix}_exp_counts.txt",
      exp_TPM             = "16_isoformswitchanalyzer_inputs/{prefix}_exp_TPM.txt",
      exp_annots          = "16_isoformswitchanalyzer_inputs/{prefix}_exp_annots.gtf",
      exp_design          = "16_isoformswitchanalyzer_inputs/{prefix}_exp_design.txt",
      exp_cpat2           = "16_isoformswitchanalyzer_inputs/{prefix}_exp_cpat2.txt",
      exp_cpat3           = "16_isoformswitchanalyzer_inputs/{prefix}_exp_cpat3.txt",
      exp_isoforms_aa     = "16_isoformswitchanalyzer_inputs/{prefix}_exp_isoforms.faa",
      exp_isoforms_nt     = "16_isoformswitchanalyzer_inputs/{prefix}_exp_isoforms-nt.fasta"
   threads:
      2
   conda:
      "envs/omics-pipelines.yaml"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("16_isoformswitchanalyzer_inputs"),
      prefix              = config["prefix"],
      refgenome_fasta     = config["refgenome_fasta"],
      isoswitch_min_count = config["isoswitch_min_count"],
      pfamscan_output     = "07_pfam/{prefix}_reference_reclocus_CDS_aa_pfamscan.pfamscan"
   log:
      out   = "logs/16_isoformswitchanalyzer_inputs_{prefix}.log"
   script:
      "scripts/get_isoformswitchanalyzer_input_files.sh"

