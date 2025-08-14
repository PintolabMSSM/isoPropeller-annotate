# ───────────────────────────────────────────────
# Rule: Get isoform count matrix
# ───────────────────────────────────────────────
rule isoform_count_matrix:
    message: "Copy the isoPropeller-collapse count matrix"
    input:
        isop_counts      = "{prefix}_exp.txt"
    output:
        iso_count_matrix = "01_isoform_counts/{prefix}_isoform-counts.txt"
    threads:
        2
    params:
        output_dir       = directory("01_isoform_counts"),
        prefix           = config["prefix"]
    log:
        "logs/01_isoform_counts/03_isoform_count_matrix_{prefix}.log"
    conda:
        SNAKEDIR + "envs/base-packages.yaml"
    shell:
        r'''
        (
           echo "## Get isoform count matrix ##"
           
           # Just copy the count file
           cat "{input.isop_counts}" > "{output.iso_count_matrix}"
        ) &> "{log}"
        '''
