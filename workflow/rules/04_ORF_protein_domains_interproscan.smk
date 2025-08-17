# ───────────────────────────────────────────────
# Local variables / imports
# ───────────────────────────────────────────────

wildcard_constraints:
    prefix = r"[^/]+",
    chunk_id = r"\d+",

def get_interproscan_chunk_ids(wildcards):
    chk_out_dir = checkpoints.interproscan_split_fasta.get(prefix=wildcards.prefix).output[0]
    with open(os.path.join(chk_out_dir, "chunk_ids.txt"), "r") as f:
        return [x for x in f.read().strip().split("\n") if x]

# ───────────────────────────────────────────────
# Checkpoint: Split protein FASTA for parallel InterProScan
# ───────────────────────────────────────────────
checkpoint interproscan_split_fasta:
    message: "Split protein FASTA into chunks for InterProScan → {wildcards.prefix}"
    input:
        fasta = f"{ANNOTATE_SUBDIRS['cds']}/{{prefix}}_reference_reclocus_CDS_aa.fa",
    output:
        directory(f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/chunks")
    params:
        seqs_per_chunk = config["seqs_per_chunk"]
    conda:
        SNAKEDIR + "envs/interproscan.yaml"
    script:
        SNAKEDIR + "scripts/split_fasta.py"


# ───────────────────────────────────────────────
# Rule: Run interproscan.sh on each protein chunk
# ───────────────────────────────────────────────
rule interproscan_per_chunk:
    message: "InterProScan: {wildcards.prefix} chunk {wildcards.chunk_id} (threads={threads})"
    input:
        faa = f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/chunks/chunk_{{chunk_id}}.fasta"
    output:
        tsv = f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.tsv",
        gff = f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gff3",
        xml = f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.xml",
    params:
        interproscan_exe = config["interproscan_path"],
        applications = config.get("interproscan_apps", "Pfam,SUPERFAMILY,MobiDBLite"),
        extra_flags  = config.get("interproscan_extra_flags", "--goterms --pathways"),
    threads: 12
    log:
        f"{INTERPROSCAN_LOG_DIR}/{{prefix}}/interproscan/{{chunk_id}}.log"
    conda:
        f"{SNAKEDIR}/envs/interproscan.yaml"
    shell:
        r'''
        (
            echo "## InterProScan chunk {wildcards.chunk_id} ({wildcards.prefix})"
            echo "Input FASTA: {input.faa}"
            echo "Using InterProScan executable: {params.interproscan_exe}"
            
            WORK_DIR=$(mktemp -d)
            
            # Use the full path to the executable from the params
            "{params.interproscan_exe}" \
                -i "{input.faa}" \
                -d "$WORK_DIR" \
                -cpu {threads} \
                -f TSV,GFF3,XML \
                -appl {params.applications} \
                {params.extra_flags}

            # Move the generated files to their final destinations
            INPUT_BASENAME=$(basename "{input.faa}")
            mv "$WORK_DIR/${{INPUT_BASENAME}}.tsv" "{output.tsv}"
            mv "$WORK_DIR/${{INPUT_BASENAME}}.gff3" "{output.gff}"
            mv "$WORK_DIR/${{INPUT_BASENAME}}.xml" "{output.xml}"

            # Sanity checks
            test -s "{output.tsv}"
            test -s "{output.gff}"
            test -s "{output.xml}"

            rm -r "$WORK_DIR"
        ) &> "{log}"
        '''

# ───────────────────────────────────────────────
# Rule: Aggregate all InterProScan chunk results
# ───────────────────────────────────────────────
rule interproscan_aggregate:
    message: "Aggregate InterProScan results → {wildcards.prefix}"
    input:
        tsv_chunks = lambda wc: expand(
            f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.tsv",
            prefix=wc.prefix, chunk_id=get_interproscan_chunk_ids(wc)
        ),
        gff_chunks = lambda wc: expand(
            f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.gff3",
            prefix=wc.prefix, chunk_id=get_interproscan_chunk_ids(wc)
        ),
        xml_chunks = lambda wc: expand(
            f"{INTERPROSCAN_OUT_DIR}/{{prefix}}/predicted/chunk_{{chunk_id}}.xml",
            prefix=wc.prefix, chunk_id=get_interproscan_chunk_ids(wc)
        ),
    output:
        tsv = f"{INTERPROSCAN_OUT_DIR}/merged/{{prefix}}.tsv",
        gff = f"{INTERPROSCAN_OUT_DIR}/merged/{{prefix}}.gff3",
        xml = f"{INTERPROSCAN_OUT_DIR}/merged/{{prefix}}.xml",
    log:
        f"{INTERPROSCAN_LOG_DIR}/aggregate/{{prefix}}.log"
    conda:
        SNAKEDIR + "envs/interproscan.yaml"
    shell:
        r'''
        (
            # echo "## Aggregating InterProScan for {wildcards.prefix}"

            # Aggregate TSV and GFF3 files (simple concatenation)
            cat {input.tsv_chunks} > "{output.tsv}"
            cat {input.gff_chunks} > "{output.gff}"

            # Aggregate XML files (requires header/footer handling)
            echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > "{output.xml}"
            echo "<protein-matches xmlns=\"http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5\">" >> "{output.xml}"
            
            # Extract protein entries from each chunk's XML file
            for xml_file in {input.xml_chunks}; do
                grep '<protein' -A 10000 "$xml_file" | grep -v '<?xml' | grep -v '<protein-matches' | grep -v '</protein-matches' >> "{output.xml}"
            done
            
            echo "</protein-matches>" >> "{output.xml}"

            # Sanity checks
            test -s "{output.tsv}"
            test -s "{output.gff}"
            test -s "{output.xml}"
        ) &> "{log}"
        '''
