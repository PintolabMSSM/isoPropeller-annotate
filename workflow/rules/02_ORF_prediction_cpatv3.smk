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

            # Use bedtools getfasta to extract isoform sequences (strand-specific, spliced)
            bedtools getfasta \
                -fi  "{input.refgenome_fasta}" \
                -bed "{output.isop_bed}" \
                -fo  "{output.isop_fasta}" \
                -s -split -name
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# Rule: Run cpatv3 with best-hit by probability
# ───────────────────────────────────────────────
rule run_cpat3_prob:
    message: "Run cpat3 with best-hit by probability"
    input:
        isop_fasta       = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
    output:
        cpat_prob_seqs   = "02_ORF_prediction/cpat_prob/{prefix}_corrected.cpatv3p18.ORF_seqs.fa",
        cpat_prob        = "02_ORF_prediction/cpat_prob/{prefix}_corrected.cpatv3p18.ORF_prob.tsv",
        cpat_prob_best   = "02_ORF_prediction/cpat_prob/{prefix}_corrected.cpatv3p18.ORF_prob.best.tsv",
        cpat_prob_no_orf = "02_ORF_prediction/cpat_prob/{prefix}_corrected.cpatv3p18.no_ORF.txt",
        cpat_prob_log    = "02_ORF_prediction/cpat_prob/{prefix}_corrected.cpatv3p18.log",
    threads: 4
    params:
        cpat_outprefix      = "02_ORF_prediction/cpat_prob/{prefix}_corrected.cpatv3p18",
        cpat_prebuilt_logit = config["cpat_prebuilt_logit"],
        cpat_prebuilt_hex   = config["cpat_prebuilt_hex"],
        minorf_cpat3        = config["minorf_cpat3"],
    resources:
        tmpdir = config["tmpdir"]
    log:
        "logs/02_ORF_prediction/cpat_prob_{prefix}.log"
    conda:
        SNAKEDIR + "envs/cpat3.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Run CPAT v3 (best-hit by PROBABILITY) ##"

            cpat.py \
               -d "{params.cpat_prebuilt_logit}" \
               -x "{params.cpat_prebuilt_hex}" \
               -g "{input.isop_fasta}" \
               --top-orf=5 \
               --min-orf={params.minorf_cpat3} \
               --best-orf=p \
               --log-file="{output.cpat_prob_log}" \
               -o "{params.cpat_outprefix}"
        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: Run cpatv3 with best-hit by length
# ───────────────────────────────────────────────
rule run_cpat3_len:
    message: "Run cpat3 with best-hit by length"
    input:
        isop_fasta       = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
    output:
        cpat_prob_seqs   = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_seqs.fa",
        cpat_prob        = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.tsv",
        cpat_prob_best   = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.best.tsv",
        cpat_prob_no_orf = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.no_ORF.txt",
        cpat_prob_log    = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.log",
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
        set -euo pipefail
        (
            echo "## Run CPAT v3 (best-hit by LENGTH) ##"

            cpat.py \
               -d "{params.cpat_prebuilt_logit}" \
               -x "{params.cpat_prebuilt_hex}" \
               -g "{input.isop_fasta}" \
               --top-orf=5 \
               --min-orf={params.minorf_cpat3} \
               --best-orf=l \
               --log-file="{output.cpat_prob_log}" \
               -o "{params.cpat_outprefix}"
        ) &> "{log}"
        '''

