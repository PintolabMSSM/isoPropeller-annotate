# ───────────────────────────────────────────────
# Rule: Generate final track files
# ───────────────────────────────────────────────
rule get_tracks:
    message: "get final track files"
    input:
        reclocus_cds_gtf      = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS.gtf",
        reclocus_cds_trns_txt = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_transcript.txt",
        reclocus_refined      = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_refined.txt",
        iso_count_matrix      = "01_isoform_counts/{prefix}_isoform-counts.txt",
        trackgroups           = "{prefix}.trackgroups",
    output:
        gtf_final    = "06_tracks/{prefix}_reference_reclocus_CDS_extra.gtf",
        gtf_stopfix  = "06_tracks/{prefix}_patched_extra_stopfix.gtf",
    params:
        prefix      = config["prefix"],
        output_dir  = "06_tracks"
    threads:
        24
    conda:
        SNAKEDIR + "envs/tracks.yaml"
    log:
        out   = "logs/06_tracks/tracks_{prefix}.log"
    script:
        SNAKEDIR + "scripts/get-tracks.sh"
