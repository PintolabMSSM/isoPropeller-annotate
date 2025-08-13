
#--------------------------------#
# ANALYZE ISOFORMS WITH PFAMSCAN #
#--------------------------------#
rule run_pfamscan:
   message: "Launch pfamscan jobs"
   input:
      reclocus_cds_aa     = "09_niap_asef/{prefix}_reference_reclocus_CDS_aa.fa"
   output:
      pfam_launched       = "07_pfam/{prefix}_jobs-submitted.check"
   threads:
      2
   conda:
      "envs/base-packages.yaml"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("07_pfam"),
      prefix              = config["prefix"],
      lsf_allocation      = config["project_allocation"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/07_pfamscan_{prefix}.log"
   script:
      "scripts/run-pfamscan.sh"

