# ───────────────────────────────────────────────
# Rule: run limma for correlations
# ───────────────────────────────────────────────
rule run_limma:
    message: "Running limma analysis"
    input:
        iso_count_matrix      = "01_isoform_counts/{prefix}_isoform-counts.txt",
        reclocus_cds_trns_txt = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_transcript.txt",
    output:
        limma_checkpoint      = "08_DEG/{prefix}_limma.checkpoint"
    params:
        prefix                = config["prefix"],
        output_dir            = "08_DEG",
        exp_design_file       = os.path.join(workflow.basedir, "../resources/exp-design_MAP-isoseq-bulk_selected-metadata.txt"),
        exp_contrasts_file    = os.path.join(workflow.basedir, "../resources/exp-contrasts_MAP-65-samples-bulk-single.txt"),
        exp_models_file       = os.path.join(workflow.basedir, "../resources/exp-models_MAP-isoseq-bulk_selected-metadata.txt"),
        exp_highlight_file    = os.path.join(workflow.basedir, "../resources/brain-tissue-specific-genes_final.txt")
    threads:
        4
    log:
        out   = "logs/08_DEG/DEG_{prefix}.log"
    conda:
        SNAKEDIR + "envs/omics-pipelines.yaml"
    script:
        SNAKEDIR + "scripts/run-limma.sh"
