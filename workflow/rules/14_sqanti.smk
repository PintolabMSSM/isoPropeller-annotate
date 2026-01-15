# ───────────────────────────────────────────────
# Run sqanti3 on the isocollapse output
# ───────────────────────────────────────────────
rule run_sqanti3:
    message: "Run sqanti3 on {wildcards.prefix}"
    input:
        isop_gtf      = "{prefix}.gtf",
        isop_fasta    = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
        genome_fasta  = config["refgenome_fasta"],
        motifs        = config["polya_motifs"],
        poly_bed      = config["polya_bed"],
        cage          = config["cage_peaks_refTSS"],
        intron        = config["intron_coverage"]
    output:
        outdir        = directory("14_sqanti/{prefix}"),
        report        = "14_sqanti/{prefix}/{prefix}_sqanti3_qc_report.pdf"
    threads: 12
    params:
        prefix        = "{prefix}_sqanti3"
    log:
        "logs/14_sqanti/{prefix}.log"
    container:
        "docker://anaconesalab/sqanti3:latest"
    shell:
        r'''
        (
            echo "## Dynamically Locating SQANTI3 Paths ##"
            
            # Locate the specific python with the sqanti3 environment
            SQANTI_PY=$(find /conda -name python | grep "/envs/sqanti3/" | head -n 1)
            
            # Locate the QC script
            SQANTI_SCRIPT=$(find /opt2 -name sqanti3_qc.py | head -n 1)
            
            echo "Using Python: $SQANTI_PY"
            echo "Using Script: $SQANTI_SCRIPT"

            $SQANTI_PY $SQANTI_SCRIPT {input.isop_fasta} {input.isop_gtf} {input.genome_fasta} \
                --polyA_motif_list "{input.motifs}" \
                --polyA_peak       "{input.poly_bed}" \
                --CAGE_peak        "{input.cage}" \
                --coverage         "{input.intron}" \
                --cpus             {threads} \
                --dir              "{output.outdir}" \
                --output           "{params.prefix}"
        
        ) &> "{log}"
        '''
