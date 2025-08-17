# ───────────────────────────────────────────────
# GET UNIQUE SPLICE CHAINS
# ───────────────────────────────────────────────
rule get_unique_splicechains:
    message: "Get unique splice chains"
    input:
        isop_bed          = "02_ORF_prediction/{prefix}_ORFpred-input.bed",
    output:
        splicechains_out  = "07_splicechains/{prefix}_splicechains.txt"
    threads:
        24
    log:
        out   = "logs/07_splicechains/splicechains_{prefix}.log"
    conda:
        SNAKEDIR + "envs/tracks.yaml"
    shell:
        r'''
        set -euo pipefail
        (
            echo "## Get unique splicechains ##"
            
            bed2splicechains.pl "{input.isop_bed}" | sort -t : -k1,1 -k2,2n > "{output.splicechains_out}"

        ) &> "{log}"
        '''
