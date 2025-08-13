#------------------#
# GET POISON EXONS #
#------------------#
rule get_nmd_poison_exons:
   message: "Identify NMD and poison exons"
   input:
      reclocus_cds_gtf    = "09_niap_asef/{prefix}_reference_reclocus_CDS.gtf"
   output:
      nmd_sj              = "10_nmd_poison_exons/{prefix}_reference_reclocus_CDS_nmd_sj.txt",
      nmd_sj_per_iso      = "10_nmd_poison_exons/{prefix}_reference_reclocus_CDS_nmd_sj_per_transcript.txt",
      nmd_sj_parsed       = "10_nmd_poison_exons/{prefix}_reference_reclocus_CDS_nmd_sj_per_transcript_parsed.txt"
   threads:
      24
   conda:
      "envs/base-packages.yaml"
   params:
      path                = shell_path,
      path_maptools       = shell_maptools,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("10_nmd_poison_exons"),
      prefix              = config["prefix"],
      refgenome_niap_gtf  = config["refgenome_niap_gtf"],
      nmdj_min_cds_len    = config["nmdj_min_cds_len"],
      nmdj_distance       = config["nmdj_distance"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/10_nmd_poison_exons_{prefix}.log"
   script:
      "scripts/get-nmd-poison-exons.sh"


