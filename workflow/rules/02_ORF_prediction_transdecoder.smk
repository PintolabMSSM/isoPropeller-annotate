# ───────────────────────────────────────────────
# Local variables / imports
# ───────────────────────────────────────────────

# Helper: resolve chunk ids produced by the checkpoint
def get_chunk_ids(wildcards):
    checkpoint_output = checkpoints.split_fasta.get(prefix=wildcards.prefix).output[0]
    with open(os.path.join(checkpoint_output, "chunk_ids.txt"), "r") as f:
        return f.read().strip().split("\n")


# ───────────────────────────────────────────────
# Rule: Create a single training model on the full dataset
# ───────────────────────────────────────────────
rule transdecoder_longorfs_global:
    input:
        fasta = "02_ORF_prediction/{prefix}_ORFpred-input.fasta"
    output:
        pep_file    = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/global_model/{{prefix}}_ORFpred-input.fasta.transdecoder_dir/longest_orfs.pep",
        model_dir   = directory(f"{TRANSDECODER_OUT_DIR}/{{prefix}}/global_model/{{prefix}}_ORFpred-input.fasta.transdecoder_dir")
    params:
        out_parent_dir = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/global_model",
        minorf         = config["minorf_transdecoder"]
    threads: 4
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/longorfs_global.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## TransDecoder.LongOrfs on full dataset: {wildcards.prefix}"
            mkdir -p "{params.out_parent_dir}"
            TransDecoder.LongOrfs \
                -m {params.minorf} \
                -S \
                -t "{input.fasta}" \
                --output_dir "{params.out_parent_dir}"

            # Explicitly check that the output file was created and is not empty
            test -s "{output.pep_file}"
        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Checkpoint: Split the GLOBAL longest_orfs.pep for parallel search
# ───────────────────────────────────────────────
checkpoint split_fasta:
    input:
        fasta = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/global_model/{{prefix}}_ORFpred-input.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        directory(f"{TRANSDECODER_OUT_DIR}/{{prefix}}/chunks")
    params:
        seqs_per_chunk=config["seqs_per_chunk"]
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    script:
        SNAKEDIR + "scripts/split_fasta.py"


# ───────────────────────────────────────────────
# Rule: Pfam domain search (runs on peptide chunks)
# ───────────────────────────────────────────────
rule hmmscan:
    input:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta"
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
            echo "## hmmscan Pfam chunk {wildcards.chunk_id} ({wildcards.prefix})"
            hmmscan \
                --domtblout "{output}" \
                --cpu {threads} \
                "{params.pfam_a_hmm}" \
                "{input}"

        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: DIAMOND per chunk (runs on peptide chunks)
# ───────────────────────────────────────────────
rule diamond_blastp:
    input:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta"
    output:
        f"{TRANSDECODER_OUT_DIR}/{{prefix}}/blastp/chunk_{{chunk_id}}.blastp"
    params:
        uniref90 = config["uniref90"]
    threads: 12
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/diamond_blastp/{{chunk_id}}.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## DIAMOND chunk {wildcards.chunk_id} ({wildcards.prefix})"
            
            TMP_OUT="{output}.tmp.$RANDOM"
            set +e
            diamond blastp \
                --block-size 5 \
                --sensitive \
                --query "{input}" \
                --db "{params.uniref90}" \
                --max-target-seqs 1 \
                --outfmt 6 \
                --evalue 1e-5 \
                --threads {threads} \
                --out "$TMP_OUT"
            EXIT_CODE="$?"
            set -e
        
        # Check the result
        # The job is a success if:
        #   1. The exit code is 0 (clean exit, even if no hits).
        #   OR
        #   2. The exit code is non-zero, BUT the output file is not empty (partial success).
        if [ "$EXIT_CODE" -eq 0 ] || [ -s "$TMP_OUT" ]; then
            echo "DIAMOND finished with exit code $EXIT_CODE. Output file is present."
            # On success, atomically move the temporary file to the final output
            mv "$TMP_OUT" "{output}"
            exit 0
        else
            echo "DIAMOND failed with exit code $EXIT_CODE and produced no output."
            # Clean up the empty/failed temp file
            rm -f "$TMP_OUT"
            exit 1
        fi

        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: Merge parallel homology search results
# ───────────────────────────────────────────────
rule merge_homology_results:
    input:
        blastp = lambda wc: expand(f"{TRANSDECODER_OUT_DIR}/{{prefix}}/blastp/chunk_{{chunk_id}}.blastp",
                                   prefix=wc.prefix, chunk_id=get_chunk_ids(wc)),
        pfam   = lambda wc: expand(f"{TRANSDECODER_OUT_DIR}/{{prefix}}/hmmscan/chunk_{{chunk_id}}.domtblout",
                                   prefix=wc.prefix, chunk_id=get_chunk_ids(wc))
    output:
        blastp = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/merged_homology.blastp",
        pfam   = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/merged_homology.pfam.domtblout"
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/merge_homology.log"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Merging homology evidence for {wildcards.prefix}"
            
            # Scalable merge for BLASTp results
            printf '%s\0' {input.blastp} | xargs -0 cat -- > "{output.blastp}"

            # Safely and scalably merge HMMER results
            : > "{output.pfam}"
            for f in {input.pfam}; do
                if grep -q "^#" "$f"; then
                    grep "^#" "$f" > "{output.pfam}"
                    break
                fi
            done
            
            # Scalably append all non-header lines from all chunk files
            printf '%s\0' {input.pfam} | xargs -0 grep -h -v "^#" -- >> "{output.pfam}" || true

        ) &> "{log}"
        '''


# ───────────────────────────────────────────────
# Rule: Predict once on the full dataset
# ───────────────────────────────────────────────
rule transdecoder_predict_final:
    input:
        fasta     = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
        blastp    = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/merged_homology.blastp",
        pfam      = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/merged_homology.pfam.domtblout",
        model_dir = directory(f"{TRANSDECODER_OUT_DIR}/{{prefix}}/global_model/{{prefix}}_ORFpred-input.fasta.transdecoder_dir")
    output:
        pep   = f"{TRANSDECODER_OUT_DIR}/merged/{{prefix}}.transdecoder.pep",
        cds   = f"{TRANSDECODER_OUT_DIR}/merged/{{prefix}}.transdecoder.cds",
        gff3  = f"{TRANSDECODER_OUT_DIR}/merged/{{prefix}}.transdecoder.gff3",
        bed   = f"{TRANSDECODER_OUT_DIR}/merged/{{prefix}}.transdecoder.bed"
    params:
        out_parent_dir = f"{TRANSDECODER_OUT_DIR}/{{prefix}}/global_model"
    threads: 4
    log:
        f"{TRANSDECODER_LOG_DIR}/{{prefix}}/predict_final.log"
    conda:
        SNAKEDIR + "envs/transdecoder.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## TransDecoder.Predict on full dataset {wildcards.prefix}"
            TransDecoder.Predict \
                -t "{input.fasta}" \
                --retain_pfam_hits "{input.pfam}" \
                --retain_blastp_hits "{input.blastp}" \
                -O "{params.out_parent_dir}"

            # Ensure destination exists, then move results
            mkdir -p "$(dirname "{output.pep}")"
            BASENAME=$(basename "{input.fasta}")

            mv "{params.out_parent_dir}/$BASENAME.transdecoder.pep"  "{output.pep}"
            mv "{params.out_parent_dir}/$BASENAME.transdecoder.cds"  "{output.cds}"
            mv "{params.out_parent_dir}/$BASENAME.transdecoder.gff3" "{output.gff3}"
            mv "{params.out_parent_dir}/$BASENAME.transdecoder.bed"  "{output.bed}"
        ) &> "{log}"
        '''
