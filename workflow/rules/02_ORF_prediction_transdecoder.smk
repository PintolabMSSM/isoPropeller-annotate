# ───────────────────────────────────────────────
# Local variables / imports
# ───────────────────────────────────────────────
import os

TRANSDECODER_OUT_DIR = "02_ORF_prediction/transdecoder"
TRANSDECODER_LOG_DIR = f"logs/{TRANSDECODER_OUT_DIR}"

# Helper: resolve chunk ids produced by the checkpoint
def get_chunk_ids(wildcards):
    checkpoint_output = checkpoints.split_fasta.get(prefix=wildcards.prefix).output[0]
    with open(os.path.join(checkpoint_output, "chunk_ids.txt"), "r") as f:
        return f.read().strip().split("\n")

# ───────────────────────────────────────────────
# Checkpoint: split fasta
# ───────────────────────────────────────────────
checkpoint split_fasta:
    input:
        fasta = "02_ORF_prediction/{prefix}_ORFpred-input.fasta"
    output:
        directory(f"{TRANSDECODER_OUT_DIR}/{{prefix}}/chunks")
    params:
        seqs_per_chunk = config["seqs_per_chunk"]
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    script:
        SNAKEDIR + "scripts/split_fasta.py"

# ───────────────────────────────────────────────
# Rule: Longest ORFs per chunk (writes under predicted/)
# ───────────────────────────────────────────────
rule transdecoder_longorfs:
    input:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta"
    output:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder_dir/longest_orfs.pep"
    params:
        out_parent_dir = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted",
        minorf = config["minorf_transdecoder"]
    threads: 2
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/longorfs/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## TransDecoder.LongOrfs chunk {wildcards.chunk_id} ({wildcards.prefix})"
            TransDecoder.LongOrfs \
                -m {params.minorf} \
                -S \
                -t "{input}" \
                --output_dir "{params.out_parent_dir}"
            test -s "{output}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# Rule: Pfam domain search (reads longest_orfs.pep from predicted/)
# ───────────────────────────────────────────────
rule hmmscan:
    input:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/hmmscan/chunk_{{chunk_id}}.domtblout"
    params:
        pfam_a_hmm = config["pfam_a_hmm"]
    threads: 4
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/hmmscan/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## hmmsearch Pfam chunk {wildcards.chunk_id} ({wildcards.prefix})"
            hmmsearch \
                -E 1e-10 \
                --domtblout "{output}" \
                --cpu {threads} \
                "{params.pfam_a_hmm}" "{input}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# Rule: DIAMOND per chunk (reads from predicted/.transdecoder_dir)
# ───────────────────────────────────────────────
rule diamond_blastp:
    input:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/blastp/chunk_{{chunk_id}}.blastp"
    params:
        uniref90 = config["uniref90"]
    threads: 4
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/diamond_blastp/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## DIAMOND chunk {wildcards.chunk_id} ({wildcards.prefix})"
            diamond blastp \
                --sensitive \
                --query "{input}" \
                --db "{params.uniref90}" \
                --max-target-seqs 1 \
                --outfmt 6 \
                --evalue 1e-5 \
                --threads {threads} \
                --out "{output}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# Rule: Predict per chunk (keeps -O predicted/)
# ───────────────────────────────────────────────
rule transdecoder_predict:
    input:
        fasta     = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta",
        domtblout = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/hmmscan/chunk_{{chunk_id}}.domtblout",
        blastp    = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/blastp/chunk_{{chunk_id}}.blastp"
    output:
        pep  = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder.pep",
        cds  = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder.cds",
        gff3 = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder.gff3",
        bed  = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder.bed"
    params:
        out_parent_dir = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted"
    threads: 1
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/predict/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## TransDecoder.Predict chunk {wildcards.chunk_id} ({wildcards.prefix})"
            TransDecoder.Predict \
                -t "{input.fasta}" \
                --retain_pfam_hits "{input.domtblout}" \
                --retain_blastp_hits "{input.blastp}" \
                -O "{params.out_parent_dir}"
            test -s "{output.pep}" && test -s "{output.cds}" && test -s "{output.gff3}" && test -s "{output.bed}"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# Rule: Aggregate merged outputs
# ───────────────────────────────────────────────
rule aggregate_results:
    input:
        lambda wc: expand(
            f"{TRANSDECODER_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.fasta.transdecoder.{{ext}}",
            prefix   = wc.prefix,
            chunk_id = get_chunk_ids(wc),
            ext      = wc.ext
        )
    output:
        f"{TRANSDECODER_OUT_DIR}/merged/{{prefix}}.transdecoder.{{ext}}"
    log:
        f"{TRANSDECODER_LOG_DIR}/aggregate/{{prefix}}.{{ext}}.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Aggregating .{wildcards.ext} for {wildcards.prefix}"
            if [[ "{wildcards.ext}" == "gff3" ]]; then
                head -n 1 "{input[0]}" > "{output}"
                grep -h -v '^#' {input} >> "{output}"
            else
                cat {input} > "{output}"
            fi
        ) &> "{log}"
        '''
