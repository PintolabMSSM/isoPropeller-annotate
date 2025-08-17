# ───────────────────────────────────────────────
# Local variables / imports
# ───────────────────────────────────────────────

# Keep wildcards tidy
wildcard_constraints:
    prefix = r"[^/]+"
    chunk_id = r"\d+"

def get_pfamscan_chunk_ids(wildcards):
    checkpoint_output = checkpoints.pfamscan_split_fasta.get(prefix=wildcards.prefix).output[0]
    with open(os.path.join(checkpoint_output, "chunk_ids.txt"), "r") as f:
        return f.read().strip().split("\n")


# ───────────────────────────────────────────────
# Checkpoint: Split protein FASTA for parallel pfamscan
# ───────────────────────────────────────────────
checkpoint pfamscan_split_fasta:
    message: "Split protein FASTA into chunks for PfamScan → {wildcards.prefix}"
    input:
        fasta = "02_ORF_prediction/{prefix}_ORFpred-input.fasta"
    output:
        directory(f"{PFAMSCAN_OUT_DIR}/{{prefix}}/chunks")
    params:
        seqs_per_chunk = config["seqs_per_chunk_pfam"]
    conda:
        SNAKEDIR + "envs/base-packages.yaml"
    script:
        SNAKEDIR + "scripts/split_fasta.py"


# ───────────────────────────────────────────────
# Rule: Run pfam_scan.pl on each protein chunk
# ───────────────────────────────────────────────
rule pfamscan_per_chunk:
    message: "PfamScan: {wildcards.prefix} chunk {wildcards.chunk_id} (threads={threads})"
    input:
        faa = f"{PFAMSCAN_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta"
    output:
        pfam_tbl = f"{PFAMSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.pfam.txt"
    params:
        pfam_db       = config["pfam_db"],                         # dir with Pfam-A.hmm (+ .h3* files)
        pfam_scan_exe = config.get("pfam_scan_exe", "pfam_scan.pl"),
        use_cut_ga    = config.get("pfam_mode", "cut_ga") == "cut_ga",
        e_dom         = config.get("pfam_e_dom", "1e-5"),
        e_seq         = config.get("pfam_e_seq", "1e-5"),
        extra_flags   = config.get("pfam_extra_flags", "-acc -clan_overlap"),
    threads: 8
    log:
        f"{PFAMSCAN_LOG_DIR}/{{prefix}}/pfamscan/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/pfamscan.yaml"
    shadow: "shallow"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## PfamScan chunk {wildcards.chunk_id} ({wildcards.prefix})"
            echo "Input FASTA: {input.faa}"
            echo "Pfam DB dir: {params.pfam_db}"
            mkdir -p "$(dirname "{output.pfam_tbl}")"

            # Preflight: ensure DB is hmmpress'd (h3* files present)
            if [ ! -f "{params.pfam_db}/Pfam-A.hmm.h3f" ]; then
              echo "ERROR: Pfam DB not hmmpress'd in {params.pfam_db} (missing .h3* files)" >&2
              exit 2
            fi

            if {params.use_cut_ga}; then
              mode_flags="-cut_ga"
            else
              mode_flags="-e_dom {params.e_dom} -e_seq {params.e_seq}"
            fi

            "{params.pfam_scan_exe}" \
              -fasta "{input.faa}" \
              -dir "{params.pfam_db}" \
              -cpu {threads} \
              $mode_flags \
              {params.extra_flags} \
              -outfile "{output.pfam_tbl}"

            test -s "{output.pfam_tbl}"
        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: Aggregate all pfamscan chunk results
# ───────────────────────────────────────────────
rule pfamscan_aggregate:
    message: "Aggregate PfamScan results → {wildcards.prefix}"
    input:
        pfam_chunks = lambda wc: expand(
            f"{PFAMSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.pfam.txt",
            prefix=wc.prefix, chunk_id=get_pfamscan_chunk_ids(wc)
        )
    output:
        pfam_tbl = f"{PFAMSCAN_OUT_DIR}/merged/{{prefix}}.pfam.txt"
    log:
        f"{PFAMSCAN_LOG_DIR}/aggregate/{{prefix}}.log"
    shell:
        r'''
        set -euo pipefail
        (
          echo "## Aggregating PfamScan for {wildcards.prefix}"
          mkdir -p "$(dirname "{output.pfam_tbl}")"

          first="{input.pfam_chunks[0]}"
          if grep -q '^#' "$first"; then
            grep '^#' "$first" > "{output.pfam_tbl}"
            grep -h -v '^#' {input.pfam_chunks} >> "{output.pfam_tbl}"
          else
            cat {input.pfam_chunks} > "{output.pfam_tbl}"
          fi

          test -s "{output.pfam_tbl}"
        ) &> "{log}"
        '''
