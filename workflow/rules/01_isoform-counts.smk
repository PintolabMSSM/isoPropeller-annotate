#--------------------------#
# GET ISOFORM COUNT MATRIX #
#--------------------------#
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
           
           # Set up output folder
           rm -rf   "{params.output_dir}"
           mkdir -p "{params.output_dir}"
           
           # Just copy the count file
           cat "{input.isop_counts}" > "{output.iso_count_matrix}"
        ) &> "{log}"
        '''
