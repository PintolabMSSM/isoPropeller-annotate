#-------------#
# RUN SQANTI3 #
#-------------#
rule sqanti3:
   message: "Run sqanti3 on the input bed file"
   input:
      isocollapse_gtf     = "{prefix}.gtf"
   output:
      sq3_classification  = "01_sqanti3/{prefix}_classification.txt",
      sq3_gtf             = "01_sqanti3/{prefix}_corrected.gtf.cds.gff",
      sq3_junctions       = "01_sqanti3/{prefix}_junctions.txt",
      sq3_params          = "01_sqanti3/{prefix}.params.txt",
      sq3_fasta           = "01_sqanti3/{prefix}_corrected.fasta",
      sq3_faa             = "01_sqanti3/{prefix}_corrected.faa"
   threads:
      24
   conda:
      "SQANTI3.env"       # Environment with sqanti required packages plus omics-toolkit
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("01_sqanti3"),
      sqanti3_qc          = expand("{sq3}/sqanti3_qc.py", sq3=config["sqanti_path"]),
      cdna_cupcake_path   = config["cdna_cupcake_path"],
      prefix              = config["prefix"],
      refgenome_fasta     = config["refgenome_fasta"],
      refgenome_gtf       = config["refgenome_gtf"],
      polya_motifs        = config["polya_motifs"],
      polya_bed           = config["polya_bed"],
      cage_peaks_allmerge = config["cage_peaks_allmerge"],
      intron_coverage     = config["intron_coverage"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/01_sqanti3_{prefix}.log",
      build = "logs/01_cDNA-cupcake_build_{prefix}.log"
   script:
      "scripts/run-sqanti3.sh"


#------------------------#
# CORRECT SQANTI3 OUTPUT #
#------------------------#
rule sqanti3_patching:
   message: "Patch sqanti3 output by collapsing novelgenes and fixing cage/polya distances"
   input:
      sq3_classification  = "01_sqanti3/{prefix}_classification.txt",
      sq3_gtf             = "01_sqanti3/{prefix}_corrected.gtf.cds.gff"
   output:
      sq3_classpatched    = "02_sqanti3_patched/{prefix}_patched_classification.txt",
      gtf_classpatched    = "02_sqanti3_patched/{prefix}_patched.gtf",
      gff_classpatched    = "02_sqanti3_patched/{prefix}_patched.gff",
      bed_classpatched    = "02_sqanti3_patched/{prefix}_patched.bed",
      stats_reclassify    = "02_sqanti3_patched/{prefix}_stats_sqanti-reclassify.txt",
      stats_novelgene     = "02_sqanti3_patched/{prefix}_stats_novelgene-overlaps-reclassify.txt"
   threads:
      2
   conda:
      "SQANTI3.env"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      path_maptools       = shell_maptools,
      output_dir          = directory("02_sqanti3_patched"),
      prefix              = config["prefix"],
      refgenome_fasta     = config["refgenome_fasta"],
      refgenome_gtf       = config["refgenome_gtf"],
      refgenome_anno      = config["refgenome_anno"],
      polya_bed           = config["polya_bed"],
      cage_peaks_allmerge = config["cage_peaks_allmerge"],
      intron_coverage     = config["intron_coverage"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/02_sqanti3_patching_{prefix}.log"
   script:
      "scripts/patch-sqanti3-output.sh"


#-------------------------------------#
# RERUN SQANTI REPORT ON PATCHED DATA #
#-------------------------------------#
rule sqanti3_patch_report:
   message: "Generate sqanti3 reports on patched data"
   input:
      sq3_classification  = "02_sqanti3_patched/{prefix}_patched_classification.txt",
      sq3_junctions       = "01_sqanti3/{prefix}_junctions.txt",
      sq3_params          = "01_sqanti3/{prefix}.params.txt"
   output:
      sq3_pdf_report      = "02_sqanti3_patched/{prefix}_patched_SQANTI3_report.pdf"
   threads:
      2
   conda:
      "SQANTI3.env"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      output_dir          = directory("02_sqanti3_patched"),
      sq3_fusion_patch    = os.path.join(workflow.basedir, "../resources/patch_sqanti3_report_fusion.patch"),
      prefix              = config["prefix"],
      sq3_path            = config["sqanti_path"]
   resources:
      tmpdir              = config["tmpdir"]
   log:
      out   = "logs/02_sqanti3_patch_report_{prefix}.log"
   script:
      "scripts/run-sqanti3-report.sh"

