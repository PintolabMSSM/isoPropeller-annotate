#------------------------------------#
# ANALYZE ISOFORMS WITH INTERPROSCAN #
#------------------------------------#
rule run_interproscan:
   message: "Launch interproscan jobs"
   input:
      reclocus_cds_aa     = "09_niap_asef/{prefix}_reference_reclocus_CDS_aa.fa"
   output:
      interpro_launched   = "06_interpro/{prefix}_jobs-submitted.check"
   threads:
      2
   conda:
      "envs/base-packages.yaml"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("06_interpro"),
      prefix              = config["prefix"],
      lsf_allocation      = config["project_allocation"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/06_interproscan_{prefix}.log"
   script:
      "scripts/run-interproscan.sh"
