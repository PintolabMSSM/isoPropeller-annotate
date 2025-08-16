#--------------------------#
# GET UNIQUE SPLICE CHAINS #
#--------------------------#
rule get_unique_splicechains:
   message: "Get unique splice chains"
   input:
      bed_classpatched    = "02_sqanti3_patched/{prefix}_patched.bed"
   output:
      splicechains_out    = "13_splicechains/{prefix}_patched_splicechains.txt"
   params:
      path                = shell_path,
      perl5lib            = shell_perl5lib,
      prefix              = config["prefix"],
      output_dir          = directory("13_splicechains")
   threads:
      24
   conda:
      "envs/tracks.yaml"
   log:
      out   = "logs/13_splicechains_{prefix}.log"
   script:
      "scripts/get-unique-splicechains.sh"

