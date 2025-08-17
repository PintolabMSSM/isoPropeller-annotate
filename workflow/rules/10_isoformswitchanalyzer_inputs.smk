# ───────────────────────────────────────────────
# Get input files needed for an isoformswitchanalyzer analysis
# ───────────────────────────────────────────────
rule get_isoformswitchanalyzer_input_files:
    message: "Gather the set of files that are needed as inputs for isoformswitchanalyzer"
    input:
        iso_count_matrix    = "01_isoform_counts/{prefix}_isoform-counts.txt",
        gtf_stopfix         = "06_tracks/{prefix}_patched_extra_stopfix.gtf",
        trackgroups         = "{prefix}.trackgroups",
        cpat_prob_best      = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.best.tsv",
        reclocus_cds_aa     = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_CDS_aa.fa",
        pfamscan_output     = f"{PFAMSCAN_OUT_DIR}/merged/{{prefix}}.pfam.txt"
    output:
        exp_counts          = "10_isoformswitchanalyzer_inputs/{prefix}_exp_counts.txt",
        exp_TPM             = "10_isoformswitchanalyzer_inputs/{prefix}_exp_TPM.txt",
        exp_annots          = "10_isoformswitchanalyzer_inputs/{prefix}_exp_annots.gtf",
        exp_design          = "10_isoformswitchanalyzer_inputs/{prefix}_exp_design.txt",
        exp_cpat2           = "10_isoformswitchanalyzer_inputs/{prefix}_exp_cpat2.txt",
        exp_cpat3           = "10_isoformswitchanalyzer_inputs/{prefix}_exp_cpat3.txt",
        exp_isoforms_aa     = "10_isoformswitchanalyzer_inputs/{prefix}_exp_isoforms.faa",
        exp_isoforms_nt     = "10_isoformswitchanalyzer_inputs/{prefix}_exp_isoforms-nt.fasta"
    threads:
        2
    params:
        output_dir          = "10_isoformswitchanalyzer_inputs",
        prefix              = config["prefix"],
        refgenome_fasta     = config["refgenome_fasta"],
        isoswitch_min_count = config["isoswitch_min_count"],
    log:
        out   = "logs/10_isoformswitchanalyzer_inputs/isoformswitchanalyzer-inputs_{prefix}.log"
    conda:
        SNAKEDIR + "envs/omics-pipelines.yaml"
    script:
        SNAKEDIR + "scripts/get_isoformswitchanalyzer_input_files.sh"
