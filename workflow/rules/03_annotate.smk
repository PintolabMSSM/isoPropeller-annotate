# ───────────────────────────────────────────────
# 1. Annotate Isoforms
# ───────────────────────────────────────────────
rule isop_annotate:
    message: "Annotating isoforms against reference GTF"
    input:
        isocollapse_gtf = "{prefix}.gtf",
        tss_bed         = "{prefix}_tss.bed",
    output:
        ref_gtf    = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference.gtf",
        ref_im_gtf = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im.gtf",
    params:
        refgenome_isop_gtf = config["refgenome_isop_gtf"],
    threads: 24
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_isop-annotate.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Running isoPropeller_annotate ##"

            isoPropeller_annotate \
                -q "{input.isocollapse_gtf}" -r "{params.refgenome_isop_gtf}" \
                -o "{output.ref_gtf}" -t {threads} -e "{input.tss_bed}"
            
            isoPropeller_annotate \
                -q "{input.isocollapse_gtf}" -r "{params.refgenome_isop_gtf}" \
                -o "{output.ref_im_gtf}" -t {threads} -m -e "{input.tss_bed}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 2. Create Summary Tables for Annotated GTFs
# ───────────────────────────────────────────────
rule isop_tabulate_reference:
    message: "Creating summary tables from annotated GTFs"
    input:
        ref_gtf    = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference.gtf",
        ref_im_gtf = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im.gtf",
    output:
        ref_gene_txt    = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_gene.txt",
        ref_trns_txt    = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_transcript.txt",
        ref_im_gene_txt = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im_gene.txt",
        ref_im_trns_txt = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im_transcript.txt",
    params:
        refgenome_fasta = config["refgenome_fasta"],
        intron_coverage = config["intron_coverage"],
        out_prefix      = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference",
        out_prefix_im   = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im",
    threads: 24
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_isop-tabulate-ref.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Running gtf2summary.pl on initial annotations ##"
            
            ATTRIBUTE_LIST=$(mktemp)
            printf '%s\n' gene_name gene_type status asm_gene_id ref_transcript_id \
                ref_gene_id ref_gene_name ref_gene_type > "$ATTRIBUTE_LIST"

            gtf2summary.pl -i "{input.ref_gtf}" -o "{params.out_prefix}" \
                -a "$ATTRIBUTE_LIST" -g "{params.refgenome_fasta}" -j "{params.intron_coverage}" -r -t {threads}

            gtf2summary.pl -i "{input.ref_im_gtf}" -o "{params.out_prefix_im}" \
                -a "$ATTRIBUTE_LIST" -g "{params.refgenome_fasta}" -j "{params.intron_coverage}" -r -t {threads}
            
            rm "$ATTRIBUTE_LIST"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 3. Calculate Fusion Gene Ratios
# ───────────────────────────────────────────────
rule isop_calculate_fusion_ratios:
    message: "Calculating fusion gene ratios for locus reconstruction"
    input:
        ref_gtf          = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference.gtf",
        ref_trns_txt     = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_transcript.txt",
        ref_im_trns_txt  = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im_transcript.txt",
        iso_count_matrix = "01_isoform_counts/{prefix}_isoform-counts.txt"
    output:
        fusion_ratio = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_fusion_gene_ratio.txt",
        fusion_mono  = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_fusion_monoexonic_gene_id.txt",
    params:
        prefix = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}", # Note: script uses this as a prefix for other outputs
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_isop-fusion-ratio.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Preparing inputs for locus reconstruction ##"
            WORK_DIR=$(mktemp -d)

            sed 's/#TranscriptID/transcript_id/' "{input.iso_count_matrix}" > "${WORK_DIR}/count.txt"
            merge2tables.pl -t1 "{input.ref_trns_txt}" -c1 0 -t2 "${WORK_DIR}/count.txt" -c2 0 -o "${WORK_DIR}/merged.txt" -s
            head -1 "${WORK_DIR}/count.txt" | sed 's/\t/\n/g' | tail -n+2 > "${WORK_DIR}/sample.txt"
            
            fusion_gene_exp_ratio.pl "{input.ref_gtf}" "${WORK_DIR}/merged.txt" "${WORK_DIR}/sample.txt" "{params.prefix}"
            fusion_gene_monoexonic_finder.pl "{input.ref_im_trns_txt}" "{output.fusion_mono}"
            
            echo "gene_id"$'\t'"monoexonic" > "${WORK_DIR}/id.txt"
            cat "{output.fusion_mono}" | sed 's/$/\t1/' >> "${WORK_DIR}/id.txt"
            
            merge2tables.pl -t1 "{params.prefix}_fusion_gene_ratio.txt" -c1 0 -t2 "${WORK_DIR}/id.txt" -c2 0 -o "${WORK_DIR}/merged.txt" -s
            sed 's/\t$/\t0/' "${WORK_DIR}/merged.txt" > "{output.fusion_ratio}"

            rm -r "$WORK_DIR"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 4. Run Locus Reconstruction
# ───────────────────────────────────────────────
rule isop_run_reclocus:
    message: "Running R-based locus reconstruction"
    input:
        fusion_ratio  = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_fusion_gene_ratio.txt",
        ref_gtf       = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference.gtf",
        ref_im_gtf    = f"{ANNOTATE_SUBDIRS['annot']}/{{prefix}}_reference_im.gtf",
        trackgroups   = "{prefix}.trackgroups",
    output:
        reclocus_id     = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reclocus_gene_id.txt",
        reclocus_gtf    = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus.gtf",
        reclocus_im_gtf = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_im_reclocus.gtf",
    params:
        prefix             = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}",
        refgenome_isop_gtf = config["refgenome_isop_gtf"],
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_isop-reclocus.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Running reconstructed_locus.R and filtering GTFs ##"

            reconstructed_locus.R "{params.prefix}" "{input.trackgroups}"
            cat {params.prefix}_reclocus_*_gene_id.txt | sort | uniq > "{output.reclocus_id}"

            gtf_reclocus.pl "{input.ref_gtf}" "{params.refgenome_isop_gtf}" "{output.reclocus_id}" "{output.reclocus_gtf}"
            gtf_reclocus.pl "{input.ref_im_gtf}" "{params.refgenome_isop_gtf}" "{output.reclocus_id}" "{output.reclocus_im_gtf}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 5. Create Summary Tables for Reconstructed Loci
# ───────────────────────────────────────────────
rule isop_tabulate_reclocus:
    message: "Creating summary tables from reconstructed locus GTFs"
    input:
        reclocus_gtf    = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus.gtf",
        reclocus_im_gtf = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_im_reclocus.gtf",
    output:
        reclocus_gene_txt      = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus_gene.txt",
        reclocus_trns_txt      = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus_transcript.txt",
        reclocus_im_gene_txt   = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_im_reclocus_gene.txt",
        reclocus_im_trns_txt   = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_im_reclocus_transcript.txt",
    params:
        refgenome_fasta = config["refgenome_fasta"],
        intron_coverage = config["intron_coverage"],
        out_prefix      = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus",
        out_prefix_im   = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_im_reclocus",
    threads: 12
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_isop-tabulate-reclocus.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Running gtf2summary.pl on reclocus annotations ##"

            ATTRIBUTE_LIST=$(mktemp)
            printf '%s\n' gene_name gene_type status asm_gene_id ref_transcript_id \
                ref_gene_id ref_gene_name ref_gene_type reclocus > "$ATTRIBUTE_LIST"
            
            TAG_LIST=$(mktemp)
            printf '%s\n' reconstructed_locus > "$TAG_LIST"

            gtf2summary.pl -i "{input.reclocus_gtf}" -o "{params.out_prefix}" -a "$ATTRIBUTE_LIST" -d "$TAG_LIST" \
                -g "{params.refgenome_fasta}" -j "{params.intron_coverage}" -r -t {threads}

            gtf2summary.pl -i "{input.reclocus_im_gtf}" -o "{params.out_prefix_im}" -a "$ATTRIBUTE_LIST" -d "$TAG_LIST" \
                -g "{params.refgenome_fasta}" -j "{params.intron_coverage}" -r -t {threads}

            rm "$ATTRIBUTE_LIST" "$TAG_LIST"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 6. Integrate ORF Predictions
# ───────────────────────────────────────────────
rule isop_integrate_orf_predictions:
    message: "Integrating ORF predictions from GMST, CPAT, and TransDecoder"
    input:
        reclocus_gtf   = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus.gtf",
        isop_fasta     = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
        cpat_leng_seqs = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_seqs.fa",
        cpat_leng      = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.tsv",
        gmst_fnn       = "02_ORF_prediction/gmst/merged/{prefix}.gmst.fnn",
        transdecoder   = "02_ORF_prediction/transdecoder/merged/{prefix}.transdecoder.cds",
    output:
        gmst_gtf         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_GMST_CDS.gtf",
        cpat_gtf         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CPAT_CDS.gtf",
        transdecoder_gtf = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_transdecoder_CDS.gtf",
    params:
        gmst_prefix         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_GMST_CDS",
        cpat_prefix         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CPAT_CDS",
        transdecoder_prefix = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_transdecoder_CDS",
        nmdj_distance       = config["nmdj_distance"],
        tis_efficiency      = os.path.join(workflow.basedir, "../resources/TIS_efficiency.txt"),
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_isop-integrate-orf.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Converting ORF predictions to GTF format ##"

            WORK_DIR=$(mktemp -d)
            fasta-reflow.pl "{input.gmst_fnn}" > "${WORK_DIR}/gmst.fnn"
            
            GMST2gtf.pl -g "{input.reclocus_gtf}" -r "${WORK_DIR}/gmst.fnn" -f "{input.isop_fasta}" \
                -o "{params.gmst_prefix}" -n "{params.nmdj_distance}" -e "{params.tis_efficiency}"

            CPAT2gtf.pl -g "{input.reclocus_gtf}" -f "{input.isop_fasta}" -p "{input.cpat_leng}" \
                -s "{input.cpat_leng_seqs}" -o "{params.cpat_prefix}" -n "{params.nmdj_distance}" -e "{params.tis_efficiency}"
            
            transdecoder2gtf.pl -g "{input.reclocus_gtf}" -f "{input.isop_fasta}" -r "{input.transdecoder}" \
                -o "{params.transdecoder_prefix}" -n "{params.nmdj_distance}" -e "{params.tis_efficiency}"
            
            rm -r "$WORK_DIR"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 7. Select Best CDS and Finalize GTF
# ───────────────────────────────────────────────
rule isop_select_best_cds:
    message: "Selecting best CDS prediction and creating final GTF"
    input:
        reclocus_gtf     = f"{ANNOTATE_SUBDIRS['reclocus']}/{{prefix}}_reference_reclocus.gtf",
        gmst_gtf         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_GMST_CDS.gtf",
        cpat_gtf         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CPAT_CDS.gtf",
        transdecoder_gtf = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_transdecoder_CDS.gtf",
    output:
        reclocus_cds_gtf = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS.gtf",
    params:
        min_prob_cpat3      = config["min_prob_cpat3"],
        gmst_prefix         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_GMST_CDS",
        cpat_prefix         = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CPAT_CDS",
        transdecoder_prefix = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_transdecoder_CDS",
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_isop-select-cds.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Selecting the best CDS prediction among all tools ##"
            WORK_DIR=$(mktemp -d)

            printf '%s\t%s\n' \
                "{params.gmst_prefix}" "GMST" \
                "{params.cpat_prefix}" "CPAT" \
                "{params.transdecoder_prefix}" "Transdecoder" > "${WORK_DIR}/prediction_list.txt"
            
            select_CDS_prediction.pl -i "${WORK_DIR}/prediction_list.txt" -o "${WORK_DIR}/temp_cds.gtf" -p "{params.min_prob_cpat3}"
            
            cat "{input.reclocus_gtf}" "${WORK_DIR}/temp_cds.gtf" | sort -k1,1 -k4,4n > "{output.reclocus_cds_gtf}"

            rm -r "$WORK_DIR"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 8. Final Tabulation and AA Sequence Generation
# ───────────────────────────────────────────────
rule isop_tabulate_final_cds:
    message: "Creating final summary tables and amino acid sequences"
    input:
        reclocus_cds_gtf = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS.gtf",
    output:
        reclocus_cds_trns_txt = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_transcript.txt",
        reclocus_cds_gene_txt = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_gene.txt",
        reclocus_cds_aa       = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_CDS_aa.fa",
    params:
        refgenome_fasta = config["refgenome_fasta"],
        intron_coverage = config["intron_coverage"],
        nmdj_distance   = config["nmdj_distance"],
        out_prefix      = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS",
        final_out_prefix = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_CDS",
    threads: 12
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_isop-tabulate-final.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Running final tabulation and generating AA sequences ##"

            ATTRIBUTE_LIST=$(mktemp)
            printf '%s\n' gene_name gene_type status asm_gene_id ref_gene_id ref_transcript_id \
                ref_gene_name ref_gene_type reclocus cds_source cds_support \
                cds_diversity tis_efficiency > "$ATTRIBUTE_LIST"

            TAG_LIST=$(mktemp)
            printf '%s\n' reconstructed_locus > "$TAG_LIST"

            gtf2summary.pl -i "{input.reclocus_cds_gtf}" -o "{params.out_prefix}" \
                -a "$ATTRIBUTE_LIST" -d "$TAG_LIST" -n "{params.nmdj_distance}" -s \
                -g "{params.refgenome_fasta}" -j "{params.intron_coverage}" -r -t {threads}

            # Generate AA fasta into its final destination folder
            gtf2fasta_CDS.pl -i "{input.reclocus_cds_gtf}" -o "{params.final_out_prefix}" \
                -g "{params.refgenome_fasta}" -t {threads}
            
            rm "$ATTRIBUTE_LIST" "$TAG_LIST"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 9. Run ASEF Analysis
# ───────────────────────────────────────────────
rule asef_analysis:
    message: "Running ASEF analysis for alternative splicing"
    input:
        reclocus_cds_gtf = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS.gtf",
        tss_bed          = "{prefix}_tss.bed",
        tts_bed          = "{prefix}_tts.bed",
    output:
        reclocus_NE_gtf  = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NE_exon.gtf",
        reclocus_NCE_gtf = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NCE_cds.gtf",
        ne_trns_txt      = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NE_transcript.txt",
    params:
        refgenome_isop_gtf = config["refgenome_isop_gtf"],
        promoter_width     = config["promoter_width"],
        intron_coverage    = config["intron_coverage"],
        prefix_FE          = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_FE",
        prefix_LE          = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_LE",
        prefix_NE          = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NE",
        prefix_NCE         = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NCE",
    threads: 12
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_asef-analysis.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Running ASEF modules ##"

            ASEF_FE.pl -q "{input.reclocus_cds_gtf}" -r "{params.refgenome_isop_gtf}" \
                -d "{params.promoter_width}" -o "{params.prefix_FE}" -e "{input.tss_bed}"
            
            ASEF_LE.pl -q "{input.reclocus_cds_gtf}" -r "{params.refgenome_isop_gtf}" \
                -o "{params.prefix_LE}" -e "{input.tts_bed}"
            
            ASEF_NE.pl -q "{input.reclocus_cds_gtf}" -r "{params.refgenome_isop_gtf}" \
                -j "{params.intron_coverage}" -i -o "{params.prefix_NE}" -t {threads}
            
            ASEF_NCE.pl -q "{input.reclocus_cds_gtf}" -p -r "{params.refgenome_isop_gtf}" \
                -o "{params.prefix_NCE}" -t {threads}
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 10. Compare Novel Exons and Novel Coding Exons
# ───────────────────────────────────────────────
rule asef_comparison:
    message: "Comparing novel exons from ASEF"
    input:
        reclocus_NCE_gtf = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NCE_cds.gtf",
        reclocus_NE_gtf  = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NE_exon.gtf",
    output:
        reclocus_NECDS = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NE_cds.txt",
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_asef-comparison.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Comparing novel coding vs. novel exons ##"
            ASEF_NCE_NE_comparison.pl \
                -c "{input.reclocus_NCE_gtf}" \
                -e "{input.reclocus_NE_gtf}" \
                -o "{output.reclocus_NECDS}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# 11. Final Classification Conversion
# ───────────────────────────────────────────────
rule isop_convert_classifications:
    message: "Merging tables and converting final classifications"
    input:
        reclocus_cds_trns_txt = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_transcript.txt",
        asef_ne_trns_txt      = f"{ANNOTATE_SUBDIRS['asef']}/{{prefix}}_reference_reclocus_CDS_NE_transcript.txt",
    output:
        reclocus_refined = f"{ANNOTATE_SUBDIRS['final']}/{{prefix}}_reference_reclocus_refined.txt",
    conda: 
        SNAKEDIR + "envs/isopropeller.yaml"
    log: 
        f"logs/{ANNOTATE_SUBDIRS['final']}/{{prefix}}_isop-convert-class.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Merging isoPropeller and ASEF tables for final classification ##"

            WORK_DIR=$(mktemp -d)

            merge2tables.pl -t1 "{input.reclocus_cds_trns_txt}" -c1 0 \
                -t2 "{input.asef_ne_trns_txt}" -c2 0 -o "${WORK_DIR}/merged.txt" -s
            
            classification_conversion.pl "${WORK_DIR}/merged.txt" "{output.reclocus_refined}"
            
            rm -r "$WORK_DIR"
        ) &> "{log}"
        '''
