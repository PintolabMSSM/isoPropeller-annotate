#----------------------------#
# RUN LIMMA FOR CORRELATIONS #
#----------------------------#
rule run_limma:
   message: "Running limma analysis"
   input:
      iso_count_matrix      = "03_isoform_counts/{prefix}_isoform-counts.txt",
      reclocus_cds_trns_txt = "09_niap_asef/{prefix}_reference_reclocus_CDS_transcript.txt"
   output:
      limma_checkpoint      = "14_DEG/{prefix}_limma.checkpoint"
   params:
      prefix                = config["prefix"],
      output_dir            = directory("14_DEG"),
      exp_design_file       = os.path.join(workflow.basedir, "../resources/exp-design_MAP-isoseq-bulk_selected-metadata.txt"),
      exp_contrasts_file    = os.path.join(workflow.basedir, "../resources/exp-contrasts_MAP-65-samples-bulk-single.txt"),
      exp_models_file       = os.path.join(workflow.basedir, "../resources/exp-models_MAP-isoseq-bulk_selected-metadata.txt"),
      exp_highlight_file    = os.path.join(workflow.basedir, "../resources/brain-tissue-specific-genes_final.txt")
   threads:
      4
   conda:
      "envs/omics-pipelines.yaml"
   log:
      out   = "logs/14_DEG_{prefix}.log"
   script:
      "scripts/run-limma.sh"
