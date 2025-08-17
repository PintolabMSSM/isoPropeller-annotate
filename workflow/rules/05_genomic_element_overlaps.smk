# ───────────────────────────────────────────────
# Rule: perform overlaps with genomic elements
# ───────────────────────────────────────────────
rule run_element_overlaps:
    message: "Process overlaps with genomic elements"
    input:
        reclocus_cds_gtf  = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS.gtf",
    output:
        combined_elements = "05_genomic_element_overlaps/{prefix}_genomic_element_overlaps.txt",
        sv_elements       = "05_genomic_element_overlaps/{prefix}_SV_overlaps.txt",
    threads:
        2
    params:
        output_dir   = directory("05_genomic_element_overlaps"),
        prefix       = config["prefix"],
        rmsk_bed     = config["refgenome_rmsk_bed"],
        rmsk_sel_bed = config["refgenome_rmsk_sel_bed"],
        ultracons    = config["refgenome_ultracons"],
        phylocsf_v31 = config["phylocsf_v31"],
        phylocsf_v35 = config["phylocsf_v35"],
        segdups_bed  = config["segdups_bed"],
        sv_ctrl      = config["sv_ctrl"],
        sv_nneu      = config["sv_nneu"]
    resources:
        tmpdir = config["tmpdir"]
    log:
        out   = "logs/05_genomic_element_overlaps/genomic_element_overlaps_{prefix}.log"
    conda:
        SNAKEDIR + "envs/base-packages.yaml"
    script:
        SNAKEDIR + "scripts/get-genomic-element-overlaps.sh"
