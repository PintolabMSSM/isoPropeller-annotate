#----------------------------------------#
# PERFORM OVERLAPS WITH GENOMIC ELEMENTS #
#----------------------------------------#
rule run_element_overlaps:
   message: "Process overlaps with genomic elements"
   input:
      sq3_gtf             = "01_sqanti3/{prefix}_corrected.gtf.cds.gff"
   output:
      combined_elements   = "08_genomic_element_overlaps/{prefix}_genomic_element_overlaps.txt",
      sv_elements         = "08_genomic_element_overlaps/{prefix}_SV_overlaps.txt"
   threads:
      2
   conda:
      "envs/base-packages.yaml"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("08_genomic_element_overlaps"),
      prefix              = config["prefix"],
      rmsk_bed            = config["refgenome_rmsk_bed"],
      rmsk_sel_bed        = config["refgenome_rmsk_sel_bed"],
      ultracons           = config["refgenome_ultracons"],
      phylocsf_v31        = config["phylocsf_v31"],
      phylocsf_v35        = config["phylocsf_v35"],
      segdups_bed         = config["segdups_bed"],
      sv_ctrl             = config["sv_ctrl"],
      sv_nneu             = config["sv_nneu"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/08_genomic_element_overlaps_{prefix}.log"
   script:
      "scripts/get-genomic-element-overlaps.sh"
