#----------------------------------------------------------------------------#
# GET ORF FILES (all, clust, clust_contained) FOR MASS-SPEC PEPTIDE ANALYSIS #
#----------------------------------------------------------------------------#
rule run_collapseORFs:
   message: "Collapse ORFs for Mass-Spec peptide mapping"
   input:
      reclocus_cds_aa     = "09_niap_asef/{prefix}_reference_reclocus_CDS_aa.fa"
   output:
      reclocus_cds_clust  = "15_collapsed-ORFs/{prefix}_reference_reclocus_CDS_aa_clust_header_generic.faa"
   threads:
      2
   conda:
      "envs/cdhit.yaml"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("15_collapsed-ORFs"),
      prefix              = config["prefix"]
   log:
      out   = "logs/15_collapseORFs_{prefix}.log"
   script:
      "scripts/get-ORFs-for-massspec.sh"
