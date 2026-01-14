# ───────────────────────────────────────────────
# Create an augmented reference file
# ───────────────────────────────────────────────
rule create_augmented_refs:
    message: "Create augmented reference consisting of isoPropeller transcripts supplemented with reference transcripts"
    input:
        gtf_stopfix      = "06_tracks/{prefix}_patched_extra_stopfix.gtf",
        reclocus_im_gtf  = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_im_reclocus.gtf",
    output:
        aug_isopropeller = "13_augmented-reference/{prefix}_reference_im_reclocus.gtf",
        aug_splicechains = "13_augmented-reference/{prefix}_reference_im_reclocus_splicechains.gtf",
    threads:
        2
    params:
        output_dir    = "13_augmented-reference",
        isof_prefix   = MERGED_ISOFORM_PREFIX,
    log:
        out   = "logs/13_augmented-reference/augmented-reference_{prefix}.log"
    conda:
        SNAKEDIR + "envs/omics-pipelines.yaml"
    shell:
        r'''
        (
            echo "## Generating augmented references ##"
            
            # For convenience we copy the isopropeller augmented reference so that it's easier to find
            cp "{input.reclocus_im_gtf}" "{output.aug_isopropeller}"
            
            # Here we additionally drop any reference isoforms with a splicechain match to isoPropeller transcripts
            gtf2gff.pl "{output.aug_isopropeller}"         > "{output.aug_isopropeller}.tmp.gff"
            gff2bed.pl "{output.aug_isopropeller}.tmp.gff" > "{output.aug_isopropeller}.tmp.bed"
            bed2splicechains.pl "{output.aug_isopropeller}.tmp.bed" > "{output.aug_isopropeller}.tmp.splicechains"
            cut -f 2 "{output.aug_isopropeller}.tmp.splicechains" \
                    | awk -v p="{params.isof_prefix}" '$0 ~ "," && $0 ~ p {{ gsub(/,/, "\\n"); print }}' \
                    | awk -v p="{params.isof_prefix}" '$0 !~ p' \
                    | sort -u \
                    > "{output.aug_isopropeller}.tmp.drop"
            gtf-filter-attributes.pl -v -a transcript_id -m "{output.aug_isopropeller}.tmp.drop" "{output.aug_isopropeller}" > "{output.aug_splicechains}"
            rm -f "{output.aug_isopropeller}.tmp."*
        
        ) &> "{log}"
        '''
