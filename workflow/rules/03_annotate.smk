
#-------------------#
# RUN isop AND ASEF #
#-------------------#
rule run_isop_asef:
   message: "Perform isop and ASEF analysis"
   input:
      isocollapse_gtf           = "{prefix}.gtf",
      tss_bed                   = "{prefix}_tss.bed",
      tts_bed                   = "{prefix}_tts.bed",
      isop_merge_counts         = "{prefix}.counts",
      trackgroups               = "{prefix}.trackgroups",
      isop_fasta                = "02_ORF_prediction/{prefix}_ORFpred-input.fasta",
      cpat_leng_seqs            = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_seqs.fa",
      cpat_leng                 = "02_ORF_prediction/cpat_leng/{prefix}_corrected.cpatv3l18.ORF_prob.tsv",
      gmst_fnn                  = "02_ORF_prediction/gmst/merged/{prefix}.gmst.fnn",
      transdecoder              = "02_ORF_prediction/transdecoder/merged/{prefix}.transdecoder.cds",
   output:
      ref_gtf                   = "03_annotate/{prefix}_reference.gtf",
      ref_im_gtf                = "03_annotate/{prefix}_reference_im.gtf",
      ref_gene_txt              = "03_annotate/{prefix}_reference_gene.txt",
      ref_trns_txt              = "03_annotate/{prefix}_reference_transcript.txt",
      ref_im_gene_txt           = "03_annotate/{prefix}_reference_im_gene.txt",
      ref_im_trns_txt           = "03_annotate/{prefix}_reference_im_transcript.txt",
      fusion_ratio              = "03_annotate/{prefix}_fusion_gene_ratio.txt",
      fusion_mono               = "03_annotate/{prefix}_fusion_monoexonic_gene_id.txt",
      reclocus_id               = "03_annotate/{prefix}_reclocus_gene_id.txt",
      reclocus_gtf              = "03_annotate/{prefix}_reference_reclocus.gtf",
      reclocus_im_gtf           = "03_annotate/{prefix}_reference_im_reclocus.gtf",
      reclocus_gene_txt         = "03_annotate/{prefix}_reference_reclocus_gene.txt",
      reclocus_trns_txt         = "03_annotate/{prefix}_reference_reclocus_transcript.txt",
      reclocus_im_gene_txt      = "03_annotate/{prefix}_reference_im_reclocus_gene.txt",
      reclocus_im_trns_txt      = "03_annotate/{prefix}_reference_im_reclocus_transcript.txt",
      reclocus_gmst_gtf         = "03_annotate/{prefix}_reference_reclocus_GMST_CDS.gtf",
      reclocus_cpat_gtf         = "03_annotate/{prefix}_reference_reclocus_CPAT_CDS.gtf",
      reclocus_transdecoder_gtf = "03_annotate/{prefix}_reference_reclocus_transdecoder_CDS.gtf",
      reclocus_cds_gtf          = "03_annotate/{prefix}_reference_reclocus_CDS.gtf",
      reclocus_cds_trns_txt     = "03_annotate/{prefix}_reference_reclocus_CDS_transcript.txt",
      reclocus_cds_gene_txt     = "03_annotate/{prefix}_reference_reclocus_CDS_gene.txt",
      reclocus_refined          = "03_annotate/{prefix}_reference_reclocus_refined.txt",
      reclocus_NE_gtf           = "03_annotate/{prefix}_reference_reclocus_CDS_NE_exon.gtf",
      reclocus_NCE_gtf          = "03_annotate/{prefix}_reference_reclocus_CDS_NCE_cds.gtf",
      reclocus_NECDS            = "03_annotate/{prefix}_reference_reclocus_CDS_NE_cds.txt",
      reclocus_cds_aa           = "03_annotate/{prefix}_reference_reclocus_CDS_aa.fa"
   threads:
      48
   conda:
      "envs/isoPropeller.yaml"
   params:
      output_dir          = directory("03_annotate"),
      prefix              = config["prefix"],
      refgenome_fasta     = config["refgenome_fasta"],
      refgenome_isop_gtf  = config["refgenome_isop_gtf"],
      intron_coverage     = config["intron_coverage"],
      promoter_width      = config["promoter_width"],
      nmdj_distance       = config["nmdj_distance"],
      tis_efficiency      = os.path.join(workflow.basedir, "../resources/TIS_efficiency.txt"),
      min_prob_cpat3      = config["min_prob_cpat3"],
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/03_annotate_{prefix}.log"
   script:
      "scripts/run-isop-asef.sh"

