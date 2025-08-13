
#----------------------------#
# ANALYZE ISOFORMS WITH CPAT #
#----------------------------#
rule run_cpat3:
   message: "Run cpat3"
   input:
      sq3_gtf             = "01_sqanti3/{prefix}_corrected.gtf.cds.gff"
   output:
      cpat_prob_seqs      = "05_cpatv3/{prefix}_corrected.cpatv3p18.ORF_seqs.fa",
      cpat_prob           = "05_cpatv3/{prefix}_corrected.cpatv3p18.ORF_prob.tsv",
      cpat_prob_best      = "05_cpatv3/{prefix}_corrected.cpatv3p18.ORF_prob.best.tsv",
      cpat_prob_no_orf    = "05_cpatv3/{prefix}_corrected.cpatv3p18.no_ORF.txt",
      cpat_prob_log       = "05_cpatv3/{prefix}_corrected.cpatv3p18.log",
      cpat_leng_seqs      = "05_cpatv3/{prefix}_corrected.cpatv3l18.ORF_seqs.fa",
      cpat_leng           = "05_cpatv3/{prefix}_corrected.cpatv3l18.ORF_prob.tsv",
      cpat_leng_best      = "05_cpatv3/{prefix}_corrected.cpatv3l18.ORF_prob.best.tsv",
      cpat_leng_no_orf    = "05_cpatv3/{prefix}_corrected.cpatv3l18.no_ORF.txt",
      cpat_leng_log       = "05_cpatv3/{prefix}_corrected.cpatv3l18.log"
   threads:
      4
   conda:
      "envs/cpat3.yaml"
   params:
      path                = shell_path,
      path_maptools       = shell_maptools, 
      perl5lib            = shell_perl5lib,
      output_dir          = directory("05_cpatv3"),
      prefix              = config["prefix"],
      cpat_prebuilt_logit = config["cpat_prebuilt_logit"],
      cpat_prebuilt_hex   = config["cpat_prebuilt_hex"],
      minorf_cpat3        = config["minorf_cpat3"],
      refgenome_fasta     = config["refgenome_fasta"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/05_cpatv3_{prefix}.log"
   script:
      "scripts/run-cpat3.sh"

#------------------#
# RUN TRANSDECOCER #
#------------------#
rule run_transdecoder:
   message: "Run transdecoder to identify ORFs"
   input:
      sq3_fasta           = "01_sqanti3/{prefix}_corrected.fasta"
   output:
      checkpoint          = "11_transdecoder/{prefix}_transdecoder.check"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      prefix              = config["prefix"],
      pfam_a_hmm          = config["pfam_a_hmm"],
      uniprot_sprot       = config["uniprot_sprot"],
      uniref90            = config["uniref90"],
      minorf_transdecoder = config["minorf_transdecoder"],
      output_dir          = directory("11_transdecoder")
   threads:
      24
   conda:
      "envs/transdecoder.yaml"
   log:
      out   = "logs/11_transdecoder_{prefix}.log"
   script:
      "scripts/run-transdecoder.sh"
