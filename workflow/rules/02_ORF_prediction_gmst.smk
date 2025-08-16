# ───────────────────────────────────────────────
# Local variables / imports
# ───────────────────────────────────────────────
import os

GMST_OUT_DIR = "02_ORF_prediction/gmst"
GMST_LOG_DIR = f"logs/{GMST_OUT_DIR}"

def get_gmst_chunk_ids(wildcards):
    checkpoint_output = checkpoints.gmst_split_fasta.get(prefix=wildcards.prefix).output[0]
    with open(os.path.join(checkpoint_output, "chunk_ids.txt"), "r") as f:
        return f.read().strip().split("\n")


# ───────────────────────────────────────────────
# Checkpoint: split fasta (GMST-dedicated)
# ───────────────────────────────────────────────
checkpoint gmst_split_fasta:
    input:
        fasta = "02_ORF_prediction/{prefix}.fasta"
    output:
        directory(f"{GMST_OUT_DIR}/{{prefix}}/chunks")
    params:
        seqs_per_chunk = config["seqs_per_chunk"]
    conda:
        SNAKEDIR + "envs/gmst.yaml"
    script:
        SNAKEDIR + "scripts/split_fasta.py"


# ───────────────────────────────────────────────
# Rule: GMST per chunk (using Snakemake `shadow`)
# ───────────────────────────────────────────────
rule gmst_per_chunk:
    input:
        fa = f"{GMST_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta"
    output:
        pep = f"{GMST_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gmst.faa",
        cds = f"{GMST_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gmst.fnn",
        gff = f"{GMST_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gmst.gff",
    params:
        out_basename = "chunk.gmst",  # safe because each job runs in its own shadow dir
        gmst_exe     = SNAKEDIR + "../bin/gmst/gmst.pl",
        gmst_bindir  = SNAKEDIR + "../bin/gmst/",
    threads: 1
    log:
        f"{GMST_LOG_DIR}/{{prefix}}/gmst/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/gmst.yaml"
    shadow: "shallow"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## GMST (SQANTI3 params + GFF) chunk {wildcards.chunk_id} ({wildcards.prefix})"
            echo "Input : {input.fa}"
            echo "Out basename: {params.out_basename}"

            export PATH="{params.gmst_bindir}:$PATH"

            # Run in the shadow dir; Snakemake will move declared outputs after success
            perl "{params.gmst_exe}" \
                -faa \
                --strand direct \
                --fnn \
                --format GFF \
                --output "{params.out_basename}" \
                "{input.fa}"

            # Rename to the basenames of declared outputs so Snakemake picks them up
            mv -f "{params.out_basename}"     "$(basename "{output.gff}")"
            mv -f "{params.out_basename}.faa" "$(basename "{output.pep}")"
            mv -f "{params.out_basename}.fnn" "$(basename "{output.cds}")"

            # Sanity checks
            test -s "$(basename "{output.gff}")"
            test -s "$(basename "{output.pep}")"
            test -s "$(basename "{output.cds}")"
        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: Aggregate merged outputs per prefix
# ───────────────────────────────────────────────
rule gmst_aggregate:
    input:
        pep_chunks = lambda wc: expand(
            f"{GMST_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gmst.faa",
            prefix=wc.prefix, chunk_id=get_gmst_chunk_ids(wc)
        ),
        cds_chunks = lambda wc: expand(
            f"{GMST_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gmst.fnn",
            prefix=wc.prefix, chunk_id=get_gmst_chunk_ids(wc)
        ),
        gff_chunks = lambda wc: expand(
            f"{GMST_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gmst.gff",
            prefix=wc.prefix, chunk_id=get_gmst_chunk_ids(wc)
        )
    output:
        pep = f"{GMST_OUT_DIR}/merged/{{prefix}}.gmst.faa",
        cds = f"{GMST_OUT_DIR}/merged/{{prefix}}.gmst.fnn",
        gff = f"{GMST_OUT_DIR}/merged/{{prefix}}.gmst.gff3",
    log:
        f"{GMST_LOG_DIR}/aggregate/{{prefix}}.log"
    shell:
        r'''
        set -euo pipefail
        (
          echo "## Aggregating GMST outputs for {wildcards.prefix}"
          mkdir -p "{GMST_OUT_DIR}/merged"

          # Concatenate FASTAs
          cat {input.pep_chunks} > "{output.pep}"
          cat {input.cds_chunks} > "{output.cds}"

          # Merge GFF: keep one header line if present; write merged as .gff3
          first="{input.gff_chunks[0]}"
          if head -n1 "$first" | grep -q '^#'; then
            head -n1 "$first" > "{output.gff}"
            for g in {input.gff_chunks}; do
              grep -h -v '^#' "$g" >> "{output.gff}"
            done
          else
            cat {input.gff_chunks} > "{output.gff}"
          fi
        ) &> "{log}"
        '''
