#----------------------------#
# GENERATE FINAL TRACK FILES #
#----------------------------#
rule get_tracks:
   message: "get final track files"
   input:
      reclocus_cds_gtf    = "09_niap_asef/{prefix}_reference_reclocus_CDS.gtf",
      reclocus_cds_trns_txt = "09_niap_asef/{prefix}_reference_reclocus_CDS_transcript.txt",
      reclocus_refined    = "09_niap_asef/{prefix}_reference_reclocus_refined.txt",
      iso_count_matrix    = "03_isoform_counts/{prefix}_isoform-counts.txt",
      trackgroups         = "{prefix}.trackgroups"
   output:
      gtf_final           = "12_tracks/{prefix}_reference_reclocus_CDS_extra.gtf",
      gtf_stopfix         = "12_tracks/{prefix}_patched_extra_stopfix.gtf"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      prefix              = config["prefix"],
      output_dir          = directory("12_tracks")
   threads:
      24
   conda:
      "envs/tracks.yaml"
   log:
      out   = "logs/12_tracks_{prefix}.log"
   script:
      "scripts/get-tracks.sh"

