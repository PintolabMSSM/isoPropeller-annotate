
#--------------------------------------#
# RUN POGO WITH PROVIDED PEPTIDE LIST  #
#--------------------------------------#
rule run_pogo:
   message: "Map a custom mass-spec peptide input file using PoGo."
   input:
      gtf_stopfix         = "12_tracks/{prefix}_patched_extra_stopfix.gtf",
      reclocus_cds_aa     = "09_niap_asef/{prefix}_reference_reclocus_CDS_aa.fa"
   output:
      exp_counts          = "17_PoGo/{prefix}_PoGo_mm1_1MM.bed"
   threads:
      2
   conda:
      os.path.join(workflow.basedir, "envs/base-packages.yaml")
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("17_PoGo"),
      prefix              = config["prefix"],
      pogo_bin            = os.path.join(workflow.basedir, "../bin/PoGo"),
      pogo_peptides       = config["pogo_peptides"]
   log:
      out   = "logs/17_pogo_{prefix}.log"
   script:
      os.path.join(workflow.basedir, "scripts/pogo-map-peptides.sh")
