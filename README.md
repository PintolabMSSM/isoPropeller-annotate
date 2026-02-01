# isoPropeller-annotate
Snakemake pipeline to annotate isoforms identified by the isoPropeller-collapse pipeline. The starting point for the pipeline is a `.gtf` isoform file, and an `isoform count matrix` with the raw read counts for each isoform (rows) in all samples (columns).

## Installation
Installation of the external software packages required by isoPropeller-collapse is largely handled by the pipeline itself, however there are a few prerequisites that need to be in place before running the pipeline. First, the pipeline requires that a conda environment named `snakemake` is present in your environment. This can be installed on the Mount Sinai 'Minerva' cluster using the following commands. 

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
conda create -c defaults -c bioconda -c conda-forge -n snakemake snakemake
```



## Running the isoPropeller-collapse pipeline

The `isoPropeller-collapse` pipeline is developed in Snakemake and uses a standardized structure that is expected by frequent Snakemake users and allows for easy deployment in modular workflows. For convenience we also provide a `run-isoPropeller-collapse` wrapper script that simplifies execution of the pipeline. The wrapper script expects three input files per analysis with the following standardized file name convention:

```
<PREFIX>.gtf     : gtf file, produced by isoPropeller-collapse
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
MAP_chr1_1.3     0        0              4              0
```

#### Starting the pipeline
We recommend organizing each isoPropeller-collapse analysis in separate folder. When ready, the Snakemake pipeline wrapper script can be used as follows:

```
run-isoPropeller-collapse -i <PREFIX>
```

#### Arguments for the run-isoPropeller-collapse wrapper script

```
   Usage: run-isoPropeller-collapse -i <PREFIX>

   Arguments:
    -i <string>
      File prefix for the .gtf and _exp_renamed.txt file
    -D 
      Run pipeline in debug mode to show the commands that will be executed.
    -help
      This help message
```



## Overview of pipeline outputs

The isoPropeller-collapse pipeline is organized as a series of tasks, each of which has their own output folder. An overview of each task and the outputs it produces is provided below.



### 01_isoform_counts

**Description:** This directory serves as the entry point for the pipeline's expression data. It contains the standardized isoform-level quantification derived from upstream processing.

**Contents:**

- **`{prefix}_isoform-counts.txt`**: A tab-delimited matrix of raw transcript counts. This file is the primary reference for all subsequent differential expression and isoform usage calculations.



### 02_ORF_prediction

**Description:** This directory contains the intermediate and final outputs for Open Reading Frame (ORF) prediction. The pipeline employs a multi-tool approach—integrating CPAT, GeneMark-ST (GMST), and TransDecoder—to identify coding sequences within the transcript isoforms. It includes coordinate conversions (GTF to BED/GFF), sequence extraction, and homology-based evidence (Pfam/UniRef) to refine ORF calls.

**Contents:**

- **`{prefix}_ORFpred-input.fasta`**: The strand-specific transcript sequences extracted from the reference genome, used as the primary input for all ORF prediction engines.
- **`{prefix}_ORFpred-input.bed` / `.gff`**: Transcripts converted into BED12 and GFF formats to facilitate coordinate-aware ORF searching.
- **`cpat_leng/`**: Contains CPAT v3 results, focusing on the "best-hit by length" logic to identify coding potential based on sequence features.
  - `{prefix}_corrected.cpatv3l18.ORF_prob.best.tsv`: The top-scoring ORF candidates per transcript.
- **`gmst/`**: Contains GeneMark-ST predictions, processed in parallel chunks for efficiency.
  - **`merged/`**: Consolidated GMST outputs including protein sequences (`.faa`), coding sequences (`.fnn`), and structural annotations (`.gff3`).
- **`transdecoder/`**: Contains the results of the TransDecoder pipeline, which integrates intrinsic sequence properties with external homology evidence.
  - **`global_model/`**: The trained Markov model based on the longest ORFs found in the dataset.
  - **`merged/`**: The final high-confidence ORF predictions (`.pep`, `.cds`, `.gff3`, `.bed`) after filtering through Pfam domain searches (via `hmmscan`) and UniRef90 protein similarity (via `diamond`).



### 03_annotation

**Description:** This directory contains the functional annotation and structural refinement of the isoforms. It uses isoPropeller to compare assembled transcripts against a reference, reconstructs genomic loci to handle fusion genes, and integrates the ORF predictions from Step 02 to finalize the coding sequence (CDS) for each transcript.

#### subdir/01_annot (Annotation)

- **`{prefix}_reference.gtf`**: The initial annotation of isoforms against the reference genome, categorizing them relative to known gene models.
- **`{prefix}_reference_transcript.txt`**: A comprehensive summary table containing gene names, types, and reference IDs for every detected isoform.

#### subdir/02_reclocus (Locus Reconstruction)

- **`{prefix}_reference_reclocus.gtf`**: Refined transcript models where genomic loci have been reconstructed to better account for fusion genes and complex overlaps.
- **`{prefix}_fusion_gene_ratio.txt`**: Quantitative metrics used to identify and characterize fusion events and monoexonic gene structures.
- **`{prefix}_reference_reclocus_transcript.txt`**: Summary statistics specifically for the reconstructed locus models, including "reconstructed_locus" tags.

#### subdir/03_cds (Coding Sequence Integration)

- **`{prefix}_reference_reclocus_CDS.gtf`**: The final structural GTF file. It integrates the best-performing ORF prediction (selected from GMST, CPAT, and TransDecoder) into the reconstructed transcript models.
- **`{prefix}_reference_reclocus_CDS_aa.fa`**: The final amino acid sequences for all predicted ORFs, with terminal stop codons removed for compatibility with downstream functional tools.
- **`{prefix}_reference_reclocus_CDS_transcript.txt`**: A detailed table including CDS support, TIS (Translation Initiation Site) efficiency, and NMD (Nonsense-Mediated Decay) status.

#### subdir/04_asef (Alternative Splicing & Exon Features)

- **`{prefix}_reference_reclocus_CDS_NE_exon.gtf`**: Identifies **Novel Exons (NE)** not found in the reference annotation.
- **`{prefix}_reference_reclocus_CDS_NCE_cds.gtf`**: Specifically highlights **Novel Coding Exons (NCE)** that alter the protein-coding potential.
- **`{prefix}_reference_reclocus_CDS_NE_cds.txt`**: A comparison file mapping novel coding segments to their respective structural exons.

#### subdir/05_final (Refined Classifications)

- **`{prefix}_reference_reclocus_refined.txt`**: The "master" classification table. It merges isoPropeller structural categories with ASEF alternative splicing data to provide a final, high-confidence call for each transcript (e.g., FSM, ISM, NIC, NNC).

#### subdir/06_poison_exon (NMD Analysis)

- **`{prefix}_nmd_sj_parsed.txt`**: Identifies "Poison Exons"—alternative exons that introduce a premature termination codon (PTC) and likely trigger Nonsense-Mediated Decay, based on the 50-55nt rule relative to downstream splice junctions.



### 04_functional_annotation

**Description:** This folder contains the results of comprehensive protein domain and motif searches. To handle the computational load of analyzing thousands of isoforms, the protein sequences are split into smaller chunks, processed in parallel using **InterProScan** and **PfamScan**, and then consolidated into final master reports.

#### subdir/interproscan/

High-level functional characterization using the InterPro database. This tool scans ORF protein sequences against multiple domain/motif databases and provides GO terms and pathway information.

- **`chunks/`**: Temporary subsets of the amino acid FASTA used for parallelization.
- **`predicted/`**: Raw output chunks in TSV, GFF3, and XML formats.
- **`merged/`**:
  - **`{prefix}.tsv`**: A tab-delimited file mapping isoforms to protein families, domains, and functional sites.
  - **`{prefix}.gff3`**: Structural representation of the protein matches.
  - **`{prefix}.xml`**: A complete, aggregated XML record of all protein matches, reconstructed with proper header/footer tags for downstream XML parsers.

#### subdir/pfamscan/

Dedicated domain annotation using the Pfam-A database via `pfam_scan.pl`. This provides a more focused search for conserved protein families using HMMER-based models.

- **`chunks/`**: Temporary subsets of the amino acid FASTA.
- **`predicted/`**: Raw output chunks containing Pfam domain hits.
- **`merged/`**:
  - **`{prefix}.pfam.txt`**: The final aggregated Pfam report. It includes significant hits and their corresponding E-values or Gathering Thresholds (GA), with the original Pfam metadata headers preserved.



### 05_genomic_element_overlaps

**Description:** This folder contains analysis of how the final reconstructed isoforms overlap with a diverse set of genomic elements. These overlaps are used to annotate isoforms with evolutionary information (PhyloCSF), repeat content (RepeatMasker), and potential impacts from Structural Variants (SVs).

**Contents:**

- **`{prefix}_genomic_element_overlaps.txt`**: A comprehensive master table summarizing intersections with several key genomic features, including:
  - **RepeatMasker (RMSK)**: Identification of transcripts overlapping with transposable elements or low-complexity repeats.
  - **Ultraconserved Elements**: Overlaps with highly conserved genomic regions.
  - **PhyloCSF**: Comparative genomics data indicating the likelihood that a region is a conserved coding sequence across species.
  - **Segmental Duplications**: Identification of transcripts residing within highly repetitive/duplicated regions of the genome.
- **`{prefix}_SV_overlaps.txt`**: A specialized report focusing on Structural Variant (SV) context. It maps isoforms against control and non-neutral SV datasets to identify transcript structures potentially altered or created by genomic rearrangements.



### 06_tracks

**Description:** This directory contains the finalized transcript models in a variety of track formats. It decorates the structural GTF with rich metadata (biotypes, structural categories, and CDS info) and generates per-sample and per-group track files (GTF and BED) for use in genome browsers like IGV or UCSC.

**Contents:**

- **`{prefix}_reference_reclocus_CDS_extra.gtf`**: The master annotation file. It includes the standard structural information plus additional custom attributes:
  - `gene_biotype`: Protein-coding, lncRNA, etc.
  - `niap_structural_category`: The high-level classification (e.g., FSM, ISM, NIC, NNC).
  - `niap_structural_subcategory`: Detailed isoform features (e.g., mono-exon, multi-exon).
  - `cds_type`: Classification of the predicted coding sequence.
  - `exon_number`: Explicit ranking of exons (1, 2, 3...) based on strand orientation.
- **`{prefix}_patched_extra_stopfix.gtf`**: A "clean" version of the master GTF where stop codons have been removed from the CDS features to ensure compatibility with certain downstream functional tools and databases.
- **`{prefix}_sample_\*.gtf.gz` / `.bed`**: Individual track files for every sample in your dataset. The BED files are specially formatted to use the isoform count (expression level) as the "score" field, allowing for visual intensity scaling in browsers.
- **`{prefix}_group_\*.gtf.gz` / `.bed`**: Aggregated track files based on the groups defined in your `trackgroups` file. These files sum the counts across all samples in a group to represent collective expression.



#### 07_splicechains

**Description:** This directory contains the coordinate-based string representations of every transcript's splicing pattern. By converting isoforms into unique splice chains, the pipeline can easily identify transcripts that share identical intron-exon structures, regardless of their transcript or gene IDs.

**Contents:**

- **`{prefix}_splicechains.txt`**: A sorted text file where each line represents a transcript's unique splicing signature.
  - **Structure:** These chains typically consist of a sequence of genomic coordinates for all splice junctions within a transcript.
  - **Usage:** This output is critical for deduplicating transcripts across different samples, comparing assemblies to reference annotations, and identifying exactly matching isoforms in complex genomic regions.



### 08_DEG

**Description:** This folder contains the results of differential expression analysis and correlation studies. By integrating the isoform count matrix with detailed metadata and experimental designs, the pipeline identifies transcripts that are significantly regulated across different conditions, tissues, or experimental contrasts.

**Contents:**

- **`{prefix}_limma.checkpoint`**: A sentinel file indicating the successful completion of the `limma` analysis suite.
- **Differential Expression Reports**: (Generated by the script) Statistical tables typically including log-fold changes, p-values, and adjusted p-values (FDR) for each transcript across the provided experimental contrasts.
- **Tissue-Specific Highlights**: Analysis refined by the `exp_highlight_file`, specifically focusing on the expression patterns of brain-tissue-specific genes within your isoform dataset.
- **Correlation Matrices**: Statistical outputs mapping the relationships between isoform expression levels and the variables defined in your experimental models.



### 09_collapsed-ORFs

**Description:** This folder contains collapsed and clustered amino acid sequences derived from the final ORF predictions. Because multiple transcript isoforms can encode identical or highly similar protein sequences, this step uses clustering (via **CD-HIT**) to group redundant proteins, providing a non-redundant search space for mass spectrometry analysis.

**Contents:**

- **`{prefix}_reference_reclocus_CDS_aa_clust_header_generic.faa`**: The primary non-redundant protein FASTA file.
  - **Clustering:** Identical protein sequences are collapsed into single representative entries.
  - **Header Formatting:** Headers are standardized to a generic format to ensure compatibility with various mass spectrometry search engines (like MaxQuant or FragPipe).
- **Mass-Spec Ready Files**: Additional subsets of ORFs, including "all" (complete set), "clust" (representative sequences), and "clust_contained" (mapping of which isoforms belong to which cluster), providing a bridge between the transcriptome and the proteome.



### 10_isoformswitchanalyzer_inputs

**Description:** This folder contains a synchronized set of expression, structural, and functional files. The pipeline filters and renames data from previous steps (counts, ORFs, CPAT scores, and Pfam domains) to ensure they meet the input requirements for analyzing isoform switches and their potential impact on protein domains using the IsoformSwitchAnalyzeR R package.

**Contents:**

- **`{prefix}_exp_counts.txt` / `_exp_TPM.txt`**: Count and abundance matrices filtered to include only those isoforms meeting the minimum count threshold for statistical reliability.
- **`{prefix}_exp_annots.gtf`**: The structural annotation (from the "stopfix" patched GTF) used to define transcript structures and exon coordinates.
- **`{prefix}_exp_design.txt`**: A design matrix mapping samples to experimental groups based on the provided trackgroups.
- **`{prefix}_exp_cpat3.txt`**: Coding potential scores formatted specifically for IsoformSwitchAnalyzeR to evaluate if switches occur between coding and non-coding isoforms.
- **`{prefix}_exp_isoforms.faa` / `-nt.fasta`**: The amino acid and nucleotide sequences of the isoforms, used by the analyzer to predict signal peptides or other sequence-based features.



### 11_pogo

**Description:** This directory facilitates the **proteogenomic characterization** of transcript isoforms. The rules within specifically format the predicted protein sequences and their corresponding genomic locations to allow for the precise mapping of experimental mass-spec peptide data back to the genome.

**Contents:**

- **`{prefix}_PoGo_input.fasta`**: A re-formatted protein FASTA file where headers are specifically structured (`>transcript:ID gene:ID`) to allow PoGo to link peptide sequences back to the correct transcript and gene models.
- **`{prefix}_PoGo_input.gtf`**: The coordinate-accurate GTF (derived from the "stopfix" patch) used by PoGo to resolve peptide locations across splice junctions.
- **`{prefix}_PoGo_mm0.txt` / `_mm1.txt`**: Copies of the original peptide list used for iterative mapping attempts.
- **`{prefix}_PoGo_mm1_1MM.bed`**: The final mapping result. By running at multiple mismatch thresholds (0 and 1), the pipeline provides a refined view of peptide support, where the resulting BED file shows the exact genomic coordinates of physical protein evidence.



### 12_plots

**Description:** This folder contains the final graphical summaries and statistical tables of the isoform analysis. By integrating structural classifications, expression counts, and biotype data, these plots provide a global view of gene and isoform diversity, allowing for rapid assessment of the pipeline's findings.

**Contents:**

- **`{prefix}_gene-isoform-read_stats.txt`**: A master statistics table summarizing the total number of genes, isoforms, and reads processed across the various classification tiers.
- **`{prefix}_genes-isoforms-reads_by_biotype.svg/pdf`**: Visual breakdowns showing the distribution of data across different genomic biotypes (e.g., protein-coding, lncRNA, antisense).
- **`{prefix}_genes-isoforms-reads_by_supercategory.svg/pdf`**: Plots categorizing transcripts into major structural groups such as Full Splice Match (FSM), Incomplete Splice Match (ISM), Novel in Catalog (NIC), and Novel Not in Catalog (NNC).
- **`{prefix}_isoform-expression-ranking.svg/pdf`**: A cumulative or ranked visualization showing the dynamic range of isoform expression, helping to distinguish dominant isoforms from low-abundance transcripts.
- **`{prefix}_isoform-stats-breakdown.svg/pdf`**: Detailed multi-panel plots that visualize various isoform-level metrics, providing a "health check" and summary of the entire transcriptome reconstruction.



### 13_augmented-reference

**Description:** This folder contains an "augmented" version of the reference genome annotation. It combines the original reference GTF with new transcripts identified by isoPropeller, while removing redundant reference isoforms that share identical splice chains isoPropeller transcripts. This results in a non-redundant, expanded reference suitable for downstream alignment or quantification.

**Contents:**

- **`{prefix}_reference_im_reclocus.gtf`**: A combined GTF file containing both reference and novel transcript models from the reconstructed locus step.
- **`{prefix}_reference_im_reclocus_splicechains.gtf`**: A refined, non-redundant version of the augmented reference.
  - **Splice Chain Filtering:** This file is generated by identifying every transcript's unique intron-exon junction sequence (splice chain).
  - **Redundancy Removal:** Any original reference isoform that is an exact structural match to a novel isoPropeller transcript is dropped, ensuring that each unique splicing pattern is represented only once in the final annotation.



### 14_sqanti

**Description:** This optional directory contains quality control plots of the isoforms using the SQANTI3 package, which characterizes transcripts through integration with orthogonal data such as **CAGE peaks** (for 5' end validation), **PolyA motifs/peaks** (for 3' end validation), and **intron coverage** (for splice junction support).

**Contents:**

- **`{prefix}_SQANTI3_report.pdf`**: A multi-page visual report providing a detailed breakdown of the transcriptome. It includes metrics on transcript length, distance to TSS/TTS, and the distribution of structural categories (FSM, ISM, NIC, NNC, etc.).
