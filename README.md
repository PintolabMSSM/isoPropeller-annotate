# sm_post-niap-merge
Snakemake pipeline to run a series of isoform analyses after performing `NIAP_merge.pl` across multiple samples. The starting point for the pipeline is a `.gtf` file produced by NIAP_merge.pl, and an `isoform count matrix` with the raw read counts for each isoform (rows) in all samples (columns).

## Installation
Installation of the external software packages required by post-niap-merge is largely handled by the pipeline itself, however there are a few prerequisites that need to be in place before running the pipeline. First, the pipeline requires that a conda environment named `snakemake` is present in your environment. This can be installed on the Mount Sinai 'Minerva' cluster using the following commands. 

```bash
module purge all
unset PYTHONPATH
unset PERL5LIB
unset R_LIBS
module load anaconda3
conda create -c defaults -c bioconda -c conda-forge -n snakemake snakemake mamba
```

For installation in other environments we recommend miniconda, which is a free minimal installer for [conda](https://docs.conda.io/en/latest/miniconda.html). After installation and initialization of miniconda the following 'conda create' command can be used to set up the isoseq pipeline.

```
conda create -c defaults -c bioconda -c conda-forge -n snakemake snakemake mamba
```

In addition to the `snakemake` conda environment there are several git repositories with additional scripts required by the pipeline. By default the pipeline expects these repositories to be available in the `opt` subdirectory of your home directory ($HOME/opt), but it is possible to specify a different location in the Snakemake config.yaml file. The list of repositories is as follows:

```
SQANTI3:          : Follow installation notes on https://github.com/ConesaLab/SQANTI3
cdna_cupcake:     : Follow installation notes on https://github.com/ConesaLab/SQANTI3
isoseq_pipeline   : git clone git@github.com:PintolabMSSM/isoseq_pipeline.git
utility           : git clone git@bitbucket.org:hvbakel/utility.git
igb_tools         : git clone git@bitbucket.org:hvbakel/igb-tools.git
ngs_tools         : git clone git@bitbucket.org:hvbakel/ngs-tools.git
minerva_queue_lsf : git clone git@bitbucket.org:hvbakel/minerva-queue-lsf.git
```
The `minerva_queue_lsf` repository provides the `submitjob` script required for the pfamscan and interpro scan tasks in the pipeline. Both these tasks partition the isoform sequences in chunks for parallel processing on the 'Minerva' compute cluster, which uses the LSF queueing system. Other environments will require a modification of the submitjob script.

## Running the post-niap-merge pipeline
The `sm_post-niap-merge` pipeline is developed in Snakemake and uses a standardized structure that is expected by frequent Snakemake users and allows for easy deployment in modular workflows. For convenience we also provide a `run-post-niap-merge` wrapper script that simplifies execution of the pipeline. The wrapper script expects three input files per analysis with the following standardized file name convention:

```
<PREFIX>.gtf     : gtf file, produced by NIAP_merge.pl
<PREFIX>_exp_renamed.txt  : raw flnc read counts for each isoform (rows) in all samples (columns)
```

#### Countlist input file
The tab-delimited `<PREFIX>_exp_reanmed.txt` file should contain the counts for all isoforms across all samples. An example counts file is shown below:

```
transcript_id     AN08161  CBL-DP19_hifi  CBL-DP20_hifi  CBL-DP21_hifi
MAP_chr1_0_1.1   0        4              69             7
MAP_chr1_0_1.2   2        0              17             0
MAP_chr1_1_1.1   0        0              0              14
MAP_chr1_1_1.2   0        9              9              0
MAP_chr1_1.3   0        0              4              0
```

#### Starting the pipeline
We recommend organizing each post-niap-merge analysis in separate folder. When ready, the Snakemake pipeline wrapper script can be used as follows:

```
run-post-niap-merge -i <PREFIX>
```

#### Arguments for the run-post-niap-merge wrapper script

```
   Usage: run-post-niap-merge -i <PREFIX>

   Arguments:
    -i <string>
      File prefix for the .gtf and _exp_renamed.txt file
    -D 
      Run pipeline in debug mode to show the commands that will be executed.
    -help
      This help message
```


## Overview of pipeline tasks and outputs
The post-niap-merge pipeline is organized as a series of tasks, each of which has their own output folder. An overview of each tasks and the outputs it produces is provided below.

### 01_sqanti3

The first stage of the pipeline is to run SQANTI3 to generate classification and quality control files for all isoforms in the provided gtf file. The 01_sqanti3 folder contains all standard output files from sqanti3, including:

* **{prefix}_classification.txt** (classification file)
* **{prefix}_corrected.gtf.cds.gff** (GTF file of aligned isoforms)
* **{prefix}_junctions.txt** (junction file)
* **{prefix}.params.txt** (summary of sqanti3 run parameters)
* **{prefix}_corrected.fasta** (isoform fasta file, after genome correction)
* **{prefix}_corrected.faa** (Translations of ORFs detected in isoforms during the sqanti3 analysis)

For a full overview of sqanti3 outputs and file formats, see the section on [Sqanti3 main outputs](https://github.com/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC#main) in the github documentation provided in the [Sqanti3](https://github.com/ConesaLab/SQANTI3/) repository.


### 02_sqanti3_patched

The purpose of this task is to process the sqanti3 output and patch a few issues that impact downstream data analysis. The following sub-tasks are performed as part of patching process:

* **CAGE distance correction**. In reviewing the distance to CAGE peaks reported in the Sqanti3 classification file we noticed that the distances are not always as expected. In this subtask we therefore use `bedtools` to calculate the distance in basepairs between the 5' end of each isoform and the nearest CAGE tag. A distance of zero indicates that the 5' end is within a CAGE tag.
* **PolyA peak distance correction**. This step is analogous to the CAGE analysis, but now using `bedtools` to calculate the distance between the 3' end of an isoform to the nearest PolyA tag.
* **Add information on isoforms that extend known gene annotations**. In this subtask each isoform is compared to reference annotations to determine whether the isoform extends beyond the boundaries of known reference annotations. This is useful to assess potential novel 5' and 3' exons. Potential gene joins (based on _any_ exon overlap on the same strand with reference transcripts) are also indicated, and used to reclassify fusion isforms annotated by sqanti as either 'fusion_known' (there are know fusion transcripts between the gene loci in the reference annotation) or 'fusion_novel' (the reference annotation does not contain fusion transcripts between the gene loci).
* **Collapse novel gene loci**. The sqanti3 classification file assigns a distinct novelGeneID to each isoform found in intergenic regions. That means that even if two intergenic isoforms overlap, each receives a different novelGeneID. The purpose of this subtask is to identify overlaps (_any_ exon overlap on the same strand) between novelGene isoforms and collapse these into new loci that are assigned a new novelGeneID for all isoforms in the loci. As part of this step, all novelGeneIDs are renumbered.
* **Annotate fusions between known and novelGene loci**. Sqanti3 only considers known genes when evaluating potential fusion isoforms because it only compares between isoforms and the reference genome annotations. In this step, potential fusion transcripts between known and novelGenes are considered and fusions between known and novel genes (any exon boundary overlap) are reclassified and are assigned a 'fusion_novel_revised_from_<original-structural-classification>' tag.
* **For isoforms that have a single best-hit transcript, assign the geneID of the best-hit transcript**. This fixes some issues where multiple associated_geneIDs were assigned to an isoform with only a single best-hit transcript.

As part of this task, several patched output files are produced:
* **{prefix}_patched_classification.txt** (An updated classification file that includes renumbered novelGene loci, updated fusion_known/fusion_novel structural categories, and extra columns with bedtools distances (bedtools_closest_cage_peak and bedtools_closest_polyA_peak), gene extensions, and fusion evidence)
* **{prefix}_patched.gtf** (Updated gtf file that has the revised novelGeneIDs in the gene_ID tag)
* **{prefix}_patched.gff** (Patched annotation track in gff format)
* **{prefix}_patched.gff** (Patched annotation track in bed format)
* **{prefix}_stats_sqanti-reclassify.txt** (Overview of how many isoforms were reclassified during the patching process)
* **{prefix}_stats_novelgene-overlaps-reclassify.txt** 

### 03_isoform_counts
This task is a stub from the equivalent NIAP_merge.pl. Since the count matrix is provided as an input to the sm_post-niap-merge pipeline we simply copy the count matrix file here.
### 04_isoform_terminal_exons_in_segdups
This pipeline task evaluates whether there are isoforms who have one, two, three, or four terminal exons mapping in a segmental duplication and where the terminal intron spans one or more other genes. These isoforms are candidates for potentially mis-mapped terminal exons due to segmental duplications, which can cause spurious gene fusions / read-through annotations. The task produces the following outputs:

* **{prefix}_mismapped-all.txt** (Listing of PBIDs for isoforms meeting the terminal exon mismapping criteria)
* **{prefix}_mismapped-all.bed** (Track file in bed format for all isoforms meeting the mismapping criteria. Can be used for evaluation in the UCSC browser)
* **{prefix}_mismaplocus_regions.txt** (Listing of unique regions with potentially mismapped isoforms to facilitate selection in the UCSC genome browser)


### 05_cpatv3
During this step, [CPAT version 3](https://github.com/liguowang/cpat) is run on all isoforms to produce a list of annotated ORFs greater than 6 amino acids (18 nt), ORF sequences and coding potential scores. The CPAT analysis is run twice using two different settings to detect the best ORF based on probability score (`cpatv3p` output file set) or length (`cpatv3l` output file set). Each output file set contains the following main outputs:

* **{prefix}.ORF_prob.best.tsv** (ORF coordinates and coding probability for the best hits)
* **.ORF_seqs.fa** (Coding sequences for predicted ORFs)
  

### 06_interpro
This task runs [interproscan](https://www.ebi.ac.uk/interpro/about/interproscan/) on all isoform nucleotide sequences (all positive reading frames) to detect domain matches in the PFAM, SMART and Panther databases. Outputs are provided in a variety of formats for downstream analysis.


### 07_pfam
This task runs [pfamscan](https://www.ebi.ac.uk/Tools/pfa/pfamscan/) on all isoform nucleotide sequences (all positive reading frames) to detect domain matches in the PFAM database. Outputs are provided in a variety of formats for downstream analysis.

### 08_genomic_element_overlaps
This task overlaps the exons or coding sequences (CDS) with various types of genomic elements and returns the number of isoform bases that overlap each element. The output is a single tab-delimited text file with isoform IDs in the first column, followed by columns for each element and the number of overlapping bases. Overlaps with the following elements are reported in the file `{prefix}_genomic_element_overlaps.txt`:

* overlap_segdups (at isoform exon- and CDS-level)
* overlap_repeatmasker (at isoform exon- and CDS-level)
* overlap_SINE-LINE-LTR (at isoform exon- and CDS-level)
* overlap_ultraconserved (at isoform exon- and CDS-level)
* overlap_phyloCSF_novel_v31 (at isoform exon- and CDS-level)
* overlap_phyloCSF_novel_v35 (at isoform exon- and CDS-level)

In addition, this task also overlaps isoforms with common structural variants (Deletions (DEL), Duplications (DUP), or other combinations (OTH), in Europeans (eur), or all ancestries at different allele frequencies (AF: 0.0001, 0.005, or 0.01). The number of isoform bases (exon level) overlapping each category of structural variants is reported in the file `{prefix}_SV_overlaps.txt`.

### 09_niap
This task compares all isoforms with all isoforms from the reference and annotates them. The `{prefix}.gtf` is comparee with the reference, which will then generates the `{prefix}_reference.gtf` file. The `{prefix}_reference.gtf` file contains the annoated information, which will be extracted and output to `{prefix}_reference.txt` containing the following columns.

* transcript_id: The transcript ID
* gene_id: The gene ID after annotation (The `|` indicates the concatenation of gene IDs forming fusion isoform in the order of 5' -> 3')
* chr: The chromosome
* strand: The strand of the isoform
* left: The left-most position of the isoform
* right: The right-most position of the isoform
* exon: The number of exons in the isoform
* length: The length of the isoform
* fraction_downstream_A: The fraction of As downstream of the 3' end in a 20 bp window
* length_downstream_A: The maximal length of As downstream of the 3' end in a 20 bp window
* gc_content: % GC of the isoform 
* asm_gene_id: The gene ID shown in `{prefix}_merged.gtf`
* gene_name: The gene name corresponding to gene_id
* gene_type: The gene type corresponding to gene_id
* ref_gene_id: The gene ID of novel gene associated gene(s) in the reference
* ref_transcript_id: The transcript ID of the most similar reference isoform
* status: The classifcation tag

An additional output file is generated with isoform-level classification information named `{prefix}_reference_refined.txt` , containing the following columns:

- transcript_id : The transcript ID
- status: The isoform status as classified by NIAP
- niap_structural_category: The isform structural category
- niap_structural_subcategory: The isoform structural subcategory

### 10_nmd_poison_exons
This task identifies isoforms that are likely targets for nonsense-mediated decay (NMD) as well as poison exons. Poison exons are naturally occurring, highly conserved alternative exons that contain a premature termination codon. Inclusion of a poison exon in a transcript targets the transcript for NMD, decreasing the amount of protein produced. 

### 11_transdecoder
This task runs a [transdecoder](https://github.com/TransDecoder/TransDecoder/wiki) analysis to identify candidate open reading frames (ORFs) within isoform sequences. In the first stage of the analysis it extracts all long ORFs per transcript above a certain threshold (defined in `config.yaml` and set to 100 amino acids by default). In the second stage [diamond blast](https://github.com/bbuchfink/diamond) and [pfamscan](https://www.ebi.ac.uk/Tools/pfa/pfamscan/) are run for each ORF to identify known domains and blast hits to [uniref90](https://www.uniprot.org/help/uniref). In the final step transdecoder predicts the likely coding regions, retaining ORFs with homology to known proteins (identified by diamond blast) or known protein domains (identified by pfamscan). Outputs are generated in bed, gff3 and pep formats.

### 12_tracks	

This task produces an extended GTF track file `{prefix}_patched_extra_stopfix.gtf` for squanti_patched isoforms (task 02). The CDS records in this GTF file have been modified to exclude the stop codon to ensure the file is compatible with `IsoswitchAnalyzeR`. The GTF file also contains the following extra attributes that are needed for e.g. a VEP analysis:

* transcript_id: The unique PBID assigned to the transcript
* gene_id: The associated gene ID(s) assigned by the NIAP pipeline. If multiple gene IDs are assigned to the same isoform they are separated by ':' (e.g. ENSG01:ENSG02)
* gene_name: The associated gene name(s) assigned by the NIAP pipelines. If multiple gene names are assigned to the same isoform they are separated by ':' (e.g. geneA:geneB)
* gene_biotype: The assigned gene biotypes(s) derived from the associated gene ID(s). If multiple biotypes are assigned to the same isoform they are assigned
* niap_structural_category: The structural category assigned by NIAP.
* niap_structural_subcategory: The structural subcategory assigned by NIAP.
* cds_type: The CDS type assigned by NIAP.

Finally, the 12_tracks folder also contains per-sample tracks that are derived from the main  `{prefix}_patched_extra_stopfix.gtf` file after restricting to isoforms present in each sample. The selection is done by selecting isoforms with a count greater than zero in the flcount matrix (task 03).

### 13_splicechains	

This task produces a file with unique splice chains derived from the patched sqanti gtf/bed files. A splicechain is a unique chain of donor/acceptor sites that can be used to compare isoforms based on splice events only, ignoring any differences that may exist in 5' and/or 3' exon lengths.
