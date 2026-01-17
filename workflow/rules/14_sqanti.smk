# ───────────────────────────────────────────────
# Run sqanti3 on the isocollapse output
# ───────────────────────────────────────────────
rule run_sqanti3:
    message: "Run sqanti3 on {wildcards.prefix}"
    input:
        isop_gtf        = "{prefix}.gtf",
        refgenome_gtf   = config["refgenome_gtf"],
        refgenome_fasta = config["refgenome_fasta"],
        motifs          = config["polya_motifs"],
        poly_bed        = config["polya_bed"],
        cage            = config["cage_peaks_refTSS"],
        intron          = config["intron_coverage"]
    output:
        report          = "14_sqanti/{prefix}/{prefix}_sqanti3_qc_report.pdf"
    threads: 12
    params:
        outdir          = "14_sqanti/{prefix}/",
        prefix          = "{prefix}_sqanti3"
    log:
        "logs/14_sqanti/{prefix}.log"
    container:
        "docker://anaconesalab/sqanti3:latest"
    shell:
        r'''
        (
            echo "## Dynamically Locating SQANTI3 Paths ##"
            
            # Create output dir
            mkdir -p {params.outdir}        
            
            # Activate the sqanti3 environment within the docker
            . /conda/miniconda3/bin/activate sqanti3
            
            # Locate the QC script
            SQANTI_SCRIPT=$(find /opt2 -name sqanti3_qc.py | head -n 1)
            
            echo "Using Python: $(which python)"
            echo "Using TD2: $(which TD2.LongOrfs)"
            echo "Using Script: $SQANTI_SCRIPT"

            python $SQANTI_SCRIPT \
                --force_id_ignore \
                --isoforms         "{input.isop_gtf}" \
                --refGTF           "{input.refgenome_gtf}" \
                --refFasta         "{input.refgenome_fasta}" \
                --polyA_motif_list "{input.motifs}" \
                --polyA_peak       "{input.poly_bed}" \
                --CAGE_peak        "{input.cage}" \
                --coverage         "{input.intron}" \
                --report            pdf \
                -t                  {threads} \
                -d                 "{params.outdir}" \
                -o                 "{params.prefix}"
        
        ) &> "{log}"
        '''
