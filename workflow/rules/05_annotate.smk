
#-------------------#
# RUN NIAP AND ASEF #
#-------------------#
rule run_niap_asef:
   message: "Perform NIAP and ASEF analysis"
   input:
      isocollapse_gtf           = "{prefix}.gtf",
      tss_bed                   = "{prefix}_tss.bed",
      tts_bed                   = "{prefix}_tts.bed",
      niap_merge_counts         = "{prefix}.counts",
      trackgroups               = "{prefix}.trackgroups",
      sq3_fasta                 = "01_sqanti3/{prefix}_corrected.fasta",
      cpat_leng                 = "05_cpatv3/{prefix}_corrected.cpatv3l18.ORF_prob.tsv",
      cpat_leng_seqs            = "05_cpatv3/{prefix}_corrected.cpatv3l18.ORF_seqs.fa",
      transdecoder_checkpoint   = "11_transdecoder/{prefix}_transdecoder.check"
   output:
      ref_gtf                   = "09_niap_asef/{prefix}_reference.gtf",
      ref_im_gtf                = "09_niap_asef/{prefix}_reference_im.gtf",
      ref_gene_txt              = "09_niap_asef/{prefix}_reference_gene.txt",
      ref_trns_txt              = "09_niap_asef/{prefix}_reference_transcript.txt",
      ref_im_gene_txt           = "09_niap_asef/{prefix}_reference_im_gene.txt",
      ref_im_trns_txt           = "09_niap_asef/{prefix}_reference_im_transcript.txt",
      fusion_ratio              = "09_niap_asef/{prefix}_fusion_gene_ratio.txt",
      fusion_mono               = "09_niap_asef/{prefix}_fusion_monoexonic_gene_id.txt",
      reclocus_id               = "09_niap_asef/{prefix}_reclocus_gene_id.txt",
      reclocus_gtf              = "09_niap_asef/{prefix}_reference_reclocus.gtf",
      reclocus_im_gtf           = "09_niap_asef/{prefix}_reference_im_reclocus.gtf",
      reclocus_gene_txt         = "09_niap_asef/{prefix}_reference_reclocus_gene.txt",
      reclocus_trns_txt         = "09_niap_asef/{prefix}_reference_reclocus_transcript.txt",
      reclocus_im_gene_txt      = "09_niap_asef/{prefix}_reference_im_reclocus_gene.txt",
      reclocus_im_trns_txt      = "09_niap_asef/{prefix}_reference_im_reclocus_transcript.txt",
      reclocus_gmst_gtf         = "09_niap_asef/{prefix}_reference_reclocus_GMST_CDS.gtf",
      reclocus_cpat_gtf         = "09_niap_asef/{prefix}_reference_reclocus_CPAT_CDS.gtf",
      reclocus_transdecoder_gtf = "09_niap_asef/{prefix}_reference_reclocus_transdecoder_CDS.gtf",
      reclocus_cds_gtf          = "09_niap_asef/{prefix}_reference_reclocus_CDS.gtf",
      reclocus_cds_trns_txt     = "09_niap_asef/{prefix}_reference_reclocus_CDS_transcript.txt",
      reclocus_cds_gene_txt     = "09_niap_asef/{prefix}_reference_reclocus_CDS_gene.txt",
      reclocus_refined          = "09_niap_asef/{prefix}_reference_reclocus_refined.txt",
      reclocus_NE_gtf           = "09_niap_asef/{prefix}_reference_reclocus_CDS_NE_exon.gtf",
      reclocus_NCE_gtf          = "09_niap_asef/{prefix}_reference_reclocus_CDS_NCE_cds.gtf",
      reclocus_NECDS            = "09_niap_asef/{prefix}_reference_reclocus_CDS_NE_cds.txt",
      reclocus_cds_aa           = "09_niap_asef/{prefix}_reference_reclocus_CDS_aa.fa"
   threads:
      48
   conda:
      "envs/base-packages.yaml"
   params:
      path                = shell_path,
      path_maptools       = shell_maptools,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("09_niap_asef"),
      prefix              = config["prefix"],
      refgenome_fasta     = config["refgenome_fasta"],
      refgenome_niap_gtf  = config["refgenome_niap_gtf"],
      intron_coverage     = config["intron_coverage"],
      promoter_width      = config["promoter_width"],
      nmdj_distance       = config["nmdj_distance"],
      gmst_fnn            = "01_sqanti3/GMST/GMST_tmp.fnn",
      tis_efficiency      = os.path.join(workflow.basedir, "../resources/TIS_efficiency.txt"),
      min_prob_cpat3      = config["min_prob_cpat3"],
      transdecoder        = "11_transdecoder/{prefix}_corrected.fasta.transdecoder.cds"
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/09_niap_asef_{prefix}.log"
   script:
      "scripts/run-niap-asef.sh"

