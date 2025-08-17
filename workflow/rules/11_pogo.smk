# ───────────────────────────────────────────────
# Run PoGo with a provided peptide list
# ───────────────────────────────────────────────
rule run_pogo:
    message: "Map a custom mass-spec peptide input file using PoGo."
    input:
        gtf_stopfix      = "06_tracks/{prefix}_patched_extra_stopfix.gtf",
        reclocus_cds_aa  = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_CDS_aa.fa",
    output:
        exp_counts       = "11_pogo/{prefix}_PoGo_mm1_1MM.bed"
    threads:
        2
    params:
        output_dir    = directory("11_pogo"),
        prefix        = config["prefix"],
        pogo_bin      = os.path.join(workflow.basedir, "../bin/PoGo"),
        pogo_peptides = config["pogo_peptides"]
    log:
        out   = "logs/11_pogo/pogo_{prefix}.log"
    conda:
        SNAKEDIR + "envs/base-packages.yaml"
    script:
        SNAKEDIR + "scripts/pogo-map-peptides.sh"
