#----------------------------------------------------------------------------#
# GET ORF FILES (all, clust, clust_contained) FOR MASS-SPEC PEPTIDE ANALYSIS #
#----------------------------------------------------------------------------#
rule run_collapseORFs:
    message: "Collapse ORFs for Mass-Spec peptide mapping"
    input:
        reclocus_cds_aa     = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_CDS_aa.fa",
    output:
        reclocus_cds_clust  = "09_collapsed-ORFs/{prefix}_reference_reclocus_CDS_aa_clust_header_generic.faa"
    threads:
        2
    params:
        output_dir  = "09_collapsed-ORFs",
        prefix      = config["prefix"]
    log:
        out   = "logs/09_collapsed-ORFs/collapseORFs_{prefix}.log"
    conda:
        SNAKEDIR + "envs/cdhit.yaml"
    script:
        SNAKEDIR + "scripts/get-ORFs-for-massspec.sh"
