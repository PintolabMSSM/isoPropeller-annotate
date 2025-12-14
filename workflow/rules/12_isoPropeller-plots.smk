# ───────────────────────────────────────────────
# Generate isoPropeller annotate plots
# ───────────────────────────────────────────────
rule isopropeller_annotate_plots:
    message: "Create isoPropeller plots with breakdowns of isoform statistics"
    input:
        file_isop_base      = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_transcript.txt",
        file_isop_refined   = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_refined.txt",
        file_isop_counts    = f"01_isoform_counts/{PREFIX}_isoform-counts.txt",
        file_isop_tgroups   = "{prefix}.trackgroups",
    output:
        txt_statstable      = "12_plots/{prefix}_gene-isoform-read_stats.txt",
        svg_isobiotype      = "12_plots/{prefix}_genes-isoforms-reads_by_biotype.svg",
        svg_isosupercat     = "12_plots/{prefix}_genes-isoforms-reads_by_supercategory.svg",
        svg_isoranking      = "12_plots/{prefix}_isoform-expression-ranking.svg",
        svg_isostats        = "12_plots/{prefix}_isoform-stats-breakdown.svg",
        pdf_isobiotype      = "12_plots/{prefix}_genes-isoforms-reads_by_biotype.pdf",
        pdf_isosupercat     = "12_plots/{prefix}_genes-isoforms-reads_by_supercategory.pdf",
        pdf_isoranking      = "12_plots/{prefix}_isoform-expression-ranking.pdf",
        pdf_isostats        = "12_plots/{prefix}_isoform-stats-breakdown.pdf",
    threads:
        2
    params:
        output_dir    = "12_plots",
        prefix        = lambda w: w.prefix,
        plot_bin      = os.path.join(workflow.basedir, "scripts/isoPropeller-isoform-stats.R"),
    log:
        out   = "logs/12_plots/plot_{prefix}.log"
    conda:
        SNAKEDIR + "envs/annotate-ggplot2.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Generating isopropeller annotate plots ##"
            
            "{params.plot_bin}" \
                -b "{input.file_isop_base}" \
                -r "{input.file_isop_refined}" \
                -c "{input.file_isop_counts}" \
                -t "{input.file_isop_tgroups}" \
                -p "{params.output_dir}/{params.prefix}"
        
        ) &> "{log}"
        '''
