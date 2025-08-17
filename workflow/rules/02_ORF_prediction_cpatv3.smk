# ───────────────────────────────────────────────
# Rule: Create fasta file with transcript sequences (strand-specific)
# ───────────────────────────────────────────────
rule get_isoform_fasta:
    message: "Create fasta file with transcript sequences"
    input:
        isop_gtf        = "{prefix}.gtf",
        refgenome_fasta = config["refgenome_fasta"],
    output:
        isop_gff   = "02_ORF_prediction/{prefix}_ORFpred-input.gff",
        isop_bed   = "02_ORF_prediction/{prefix}_ORFpred-input.bed",
        isop_fasta = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
    threads: 2
    log:
        "logs/02_ORF_prediction/{prefix}_isoform-fasta.log"
    conda:
        SNAKEDIR + "envs/base-packages.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Create isoform fasta file with bedtools ##"

            # Get bed12 formatted isoform annotations
            gtf2gff.pl "{input.isop_gtf}"  > "{output.isop_gff}"
            gff2bed.pl "{output.isop_gff}" > "{output.isop_bed}"

            # Run bedtools and pipe its output to sed for header cleaning
            bedtools getfasta \
                -fi  "{input.refgenome_fasta}" \
                -bed "{output.isop_bed}" \
                -s -split \
                -nameOnly \
            | sed '/^>/ s/\s*(.*)//' > "{output.isop_fasta}"

        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: Run cpatv3 with best-hit by length (using BED input)
# ───────────────────────────────────────────────
rule run_cpat3_len:
    message: "Run cpat3 with best-hit by length using BED and genome"
    input:
        # Changed from transcript FASTA to BED file and genome FASTA
        isop_bed        = "02_ORF_prediction/{prefix}_ORFpred-input.bed",
        refgenome_fasta = config["refgenome_fasta"],
    output:
        cpat_leng_seqs   = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_seqs.fa",
        cpat_leng        = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.tsv",
        cpat_leng_best   = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.best.tsv",
        cpat_leng_no_orf = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.no_ORF.txt",
        cpat_leng_log    = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.log",
    threads: 4
    params:
        cpat_outprefix      = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18",
        cpat_prebuilt_logit = config["cpat_prebuilt_logit"],
        cpat_prebuilt_hex   = config["cpat_prebuilt_hex"],
        minorf_cpat3        = config["minorf_cpat3"],
    resources:
        tmpdir = config["tmpdir"]
    log:
        "logs/02_ORF_prediction/cpat_leng_{prefix}.log"
    conda:
        SNAKEDIR + "envs/cpat3.yaml"
    shell:
        r'''
        (
            echo "## Run CPAT v3 (best-hit by LENGTH) with BED input ##"

            # Find the correct CPAT executable (cpat or cpat.py)
            CPAT_CMD=$(command -v cpat || command -v cpat.py)
            if [ -z "$CPAT_CMD" ]; then
                echo "Error: Neither 'cpat' nor 'cpat.py' found." >&2; exit 1
            fi
            echo "Using CPAT command: $CPAT_CMD"

            # Run CPAT using the reference genome (-g) and a BED file of transcripts (-r)
            # This avoids the ID capitalization and suffix issues.
            "$CPAT_CMD" \
                -d "{params.cpat_prebuilt_logit}" \
                -x "{params.cpat_prebuilt_hex}" \
                -r "{input.refgenome_fasta}" \
                -g "{input.isop_bed}" \
                --top-orf=5 \
                --min-orf={params.minorf_cpat3} \
                --best-orf=l \
                --log-file="{output.cpat_leng_log}" \
                -o "{params.cpat_outprefix}"
        ) &> "{log}"
        '''
