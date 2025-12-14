#!/usr/bin/env Rscript

# 14.12.2025 13:04:43 EST

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.
library(getopt)
args = matrix(c('file_isop_base'     , 'b', 1, "character", "Path to isoPropeller _reference_reclocus_CDS_transcript.txt file",
                'file_isop_refined'  , 'r', 1, "character", "Path to isoPropeller _reference_reclocus_refined.txt file",
                'file_isop_counts'   , 'c', 1, "character", "Path to isoPropeller _exp.txt count matrix file",
                'file_isop_tgroups'  , 't', 1, "character", "Path to isoPropeller .trackgroups file",
                'out_prefix'         , 'p', 2, "character", "Name of output file",
                'help'               , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$out_prefix) ) { opt$out_prefix = "isoPropeller_plots" }

# Help message
if ( !is.null(opt$help) || is.null(opt$file_isop_base) || is.null(opt$file_isop_refined) || is.null(opt$file_isop_counts) || is.null(opt$file_isop_tgroups)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############

library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(scales)
library(ggbreak)
library(aplot)


#############
# FUNCTIONS #
#############

# Clean up the SVG format for inkscape by removing 'textlength' and 'lengthAdjust' attributes
clean_svg_file <- function( svg_file="" ){
    if (grepl("\\.svg$", svg_file, ignore.case = TRUE)){
        svg_content <- readLines(svg_file)
        svg_content <- gsub(' textLength=["\'][^"\']+["\']', '', svg_content)
        svg_content <- gsub(' lengthAdjust=["\'][^"\']+["\']', '', svg_content)
        writeLines(svg_content, svg_file)
    }
}


########
# MAIN #
########

# Define names of colums containing required data
col_base_transcript_id = "transcript_id"
col_base_gene_id       = "gene_id"
col_base_gene_name     = "gene_name"
col_base_gene_type     = "gene_type"
col_base_length        = "length"
col_base_exon_count    = "exon"
col_base_cds_type      = "cds_type"
col_base_reclocus      = "reconstructed_locus"
col_rfnd_transcript_id = "transcript_id"
col_rfnd_status        = "status"
col_rfnd_category      = "isoPropeller_structural_category"
col_rfnd_subcategory   = "isoPropeller_structural_subcategory"



##-------------------
## LOAD INPUTS
##-------------------

## Load metadata
isop_base    = read.table(opt$file_isop_base,    header=T, sep="\t", check.names=F)
isop_refined = read.table(opt$file_isop_refined, header=T, sep="\t", check.names=F)
isop_base    = isop_base[, c(col_base_transcript_id, col_base_gene_id, col_base_gene_name, col_base_gene_type, 
                             col_base_length, col_base_exon_count, col_base_cds_type, col_base_reclocus) ]
isop_refined = isop_refined[, c(col_rfnd_transcript_id, col_rfnd_status, col_rfnd_category, col_rfnd_subcategory)]
d            = merge(isop_base, isop_refined)
colnames(d)  = c("transcript_id", "gene_id", "gene_name", "gene_type", "length", "exon_count", "cds_type", "reclocus", "status", "category", "subcategory")

## Add a new supercategory column
recode_subcategory_to_supercategory <- function(x) {
    recode(x, 'known'                     = 'Annotated (Known, FSM)',
              'partial'                   = 'Annotated (Known, ISM)',
              'partial_3p_extended'       = 'Annotated (Known, ISM)',
              'partial_5p_extended'       = 'Annotated (Known, ISM)',
              'partial_extended'          = 'Annotated (Known, ISM)',
              'fusion_known'              = 'Fusion/multigene',
              'fusion_mixed'              = 'Fusion/multigene',
              'fusion_novel'              = 'Fusion/multigene',
              'antisense'                 = 'Novel gene, antisense',
              'divergent'                 = 'Novel gene, divergent',
              'intergenic_no_overlap'     = 'Novel gene, intergenic',
              'intergenic_non_pc_overlap' = 'Novel gene, sense-overlap',
              'intronic'                  = 'Novel gene, sense-overlap',
              'ncRNA_host_gene'           = 'Novel gene, sense-overlap',
              'overlapping'               = 'Novel gene, sense-overlap',
              'monoexonic'                = 'Novel, monoexonic',
              'novel_in_catalog'          = 'Novel, in catalog',
              'novel_not_in_catalog'      = 'Novel, not in catalog')
}
d$supercategory = recode_subcategory_to_supercategory(d$subcategory)
d$supercategory = factor(d$supercategory, 
                         levels=c('Annotated (Known, FSM)', 'Annotated (Known, ISM)', 'Novel, in catalog', 
                                  'Novel, not in catalog', 'Novel, monoexonic', 'Fusion/multigene', 
                                  'Novel gene, antisense', 'Novel gene, divergent', 'Novel gene, intergenic', 
                                  'Novel gene, sense-overlap') )

# Process gene biotypes into a new unified_type column
d$gene_type <- sapply(strsplit(d$gene_type, "\\|"), function(parts) {
  unique_parts <- unique(parts)
  if (length(unique_parts) > 1) {
    "Multiple"
  } else {
    unique_parts
  }
})
d$gene_type[ d$gene_type == "na"]       = "Unknown"
d$gene_type[ d$gene_type == "NA"]       = "Unknown"
d$gene_type[ is.na(d$gene_type) ]       = "Unknown"
d$gene_type[ d$gene_type == "artifact"] = "Unknown"

recode_gene_type_unified_type <- function(x) {
    recode(x, 'protein_coding'                     = 'Protein-coding',
              'lncRNA'                             = 'lncRNA',
              'Unknown'                            = 'Novel',
              'Multiple'                           = 'Multiple biotypes',
              'processed_pseudogene'               = 'Processed pseudogene',
              'unprocessed_pseudogene'             = 'Unprocessed pseudogene',
              'transcribed_processed_pseudogene'   = 'TR-PR pseudogene',
              'transcribed_unprocessed_pseudogene' = 'TR-UP pseudogene',
              'transcribed_unitary_pseudogene'     = 'TR-UN pseudogene',
              'translated_unprocessed_pseudogene'  = 'TL-UP pseudogene',
              'misc_RNA'                           = 'Misc RNA',
              'TEC'                                = 'TEC',
              'miRNA'                              = 'miRNA',
              'snoRNA'                             = 'snoRNA',
              'unitary_pseudogene'                 = 'Unitary pseudogene',
              'TR_V_pseudogene'                    = 'TR_V_pseudogene',
              'TR_C_gene'                          = 'TR_C_gene',
              'TR_V_gene'                          = 'TR_V_gene',
              'IG_C_gene'                          = 'IG_C_gene',
              'IG_V_gene'                          = 'IG_V_gene' )
}
d$unified_type = recode_gene_type_unified_type(d$gene_type)
d$unified_type = factor(d$unified_type, 
                         levels=c('Protein-coding','lncRNA','Novel','Misc RNA','miRNA','snoRNA','TEC',
                                  'Processed pseudogene','Unprocessed pseudogene', 'Unitary pseudogene',
                                  'TR-UN pseudogene','TR-PR pseudogene','TR-UP pseudogene','TL-UP pseudogene',
                                  'TR_V_pseudogene','TR_C_gene','TR_V_gene',
                                  'IG_C_gene','IG_V_gene','Multiple biotypes') )

# Read count files
counts      = read.table(opt$file_isop_counts, header=T, sep="\t", check.names=F, row.names = 1, comment.char = "")
trackgroups = read.table(opt$file_isop_tgroups, header=T, sep="\t", check.names=F)
tglist      = c( "all", unique( trackgroups[,2]) )


# Get a summary of read counts across all samples and per trackgroup
l.counts = list()
for ( tg in tglist ){
    if (tg == "all"){
        selected = colnames(counts)
    } else{
        selected = trackgroups[ trackgroups[,2] == tg, 1] 
    }
    l.counts[[ paste0( "readcount_", tg) ]] = rowSums( counts[, selected, drop=FALSE ] )
    l.counts[[ paste0( "isocount_", tg) ]]  = as.numeric( l.counts[[ paste0( "readcount_", tg) ]] > 0 )
}
d.counts = data.frame(l.counts)
colnames(d.counts) = names(l.counts)
d.counts$transcript_id = rownames(d.counts)


# Format base ggplot objects
d = merge(d, d.counts, by="transcript_id")
d = as.data.table(d)
head(d)



##---------------------------------------------------------
## ISOFORM AND READ BREAKDOWNS ACROSS STRUCTURAL CATEGORIES
##---------------------------------------------------------

pd_iso  = d[, .(genecount = uniqueN(gene_name), isocount = sum(isocount_all), readcount = sum(readcount_all)), by = supercategory]

# For horizontal bars, reverse the factors
pd_iso$supercategory = factor(pd_iso$supercategory, levels = rev(levels(pd_iso$supercategory)))

# For horizontal bars, reverse the color palette (select the default ggplot colors from scales::hue_pal
default_colors  = scales::hue_pal()( length( levels(pd_iso$supercategory) ) )
reversed_colors = rev(default_colors)

# Calculate percentages to show next to the bars
short_label <- scales::label_number(scale_cut = cut_short_scale(), accuracy = 0.1)
pd_iso <- pd_iso %>%
    mutate(
        genepct   = genecount / sum(genecount),
        genelabel = paste0(short_label(genecount), " (", scales::percent(genepct, accuracy = 1), ")"),
        isopct    = isocount / sum(isocount),
        isolabel  = paste0(short_label(isocount),  " (", scales::percent(isopct, accuracy = 1),  ")"),
        readpct   = readcount / sum(readcount),
        readlabel = paste0(short_label(readcount), " (", scales::percent(readpct, accuracy = 1),  ")")
  )

# General plot elements
bar_width  = 0.8
text_size  = 24
label_size = 5

# Plot for gene count
plot_gcount = ggplot(pd_iso, aes(y = supercategory, x = genecount, fill = supercategory)) +
    geom_col( width = bar_width ) +
    geom_text( aes(label = genelabel), hjust = -0.1, size=label_size) +
    theme_minimal() +
    theme( text        = element_text(size=text_size),
           axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                     margin = margin(t = 10, b = 10, unit = "pt")),
           axis.text.y = element_text(hjust = 0),
           legend.position = "none") +
    labs(title = "Genes", x = "# Genes", y = "") +
    scale_fill_manual(values = reversed_colors) +
    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale(), accuracy = 1),
                       expand = expansion(mult = c(0, 0.5)))

# Plot for isoform count
plot_icount = ggplot(pd_iso, aes(y = supercategory, x = isocount, fill = supercategory)) +
    geom_col( width=bar_width ) +
    geom_text( aes(label = isolabel), hjust = -0.1, size=label_size) +
    theme_minimal() +
    theme( text            = element_text(size=text_size),
           legend.position = "none",
           axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                     margin = margin(t = 10, b = 10, unit = "pt")),
           axis.text.y     = element_blank()) +
    labs(title = "Isoforms", x = "# Isoforms", y = "") +
    scale_fill_manual(values = reversed_colors) +
    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale(), accuracy = 1),
                       expand = expansion(mult = c(0, 0.55)))

# Plot for read count
plot_rcount = ggplot(pd_iso, aes(y = supercategory, x = readcount, fill = supercategory)) +
    geom_col( width = bar_width ) +
    geom_text( aes(label = readlabel), hjust = -0.1, size=label_size) +
    theme_minimal() +
    theme( text            = element_text(size=text_size),
           legend.position = "none",
           
           # --- FIX FOR TOP AXIS ---
           # Explicitly hide top axis text and ticks
           axis.text.x.top  = element_blank(),
           axis.ticks.x.top = element_blank(),
           axis.line.x.top  = element_blank(),
           # ------------------------

           axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                    margin = margin(t = 10, b = 10, unit = "pt")),
           axis.text.y     = element_blank()) +
    labs(title = "Reads", x = "# Reads", y = "", fill="Category") +
    scale_fill_manual(values = reversed_colors) +
    
    # --- SINGLE BREAK ---
    # Cut at 1M, resume at 9M
    scale_x_break(c(1000000, 9000000), scales = 2) + 

    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale(), accuracy = 1),
                       expand = expansion(mult = c(0, 0.7)))


options(repr.plot.width=20, repr.plot.height=6)
# Combine the plots using aplot instead of patchwork (|)
# widths = c(Genes, Isoforms, Reads)
# We give "Genes" 1.5x the width of the others to account for the text labels on the left.
final_combined_plot <- plot_list(
    plot_gcount, 
    plot_icount, 
    plot_rcount, 
    ncol = 3, 
    widths = c(1.6, 1, 1) 
)

# Save using ggplot2::ggsave (this handles the aplot object correctly)
ggsave(filename = paste0(opt$out_prefix, "_genes-isoforms-reads_by_supercategory.pdf"), 
       plot = final_combined_plot, 
       width = 17, height = 6)
ggsave(filename = paste0(opt$out_prefix, "_genes-isoforms-reads_by_supercategory.svg"), 
       plot = final_combined_plot, 
       width = 17, height = 6)
clean_svg_file(paste0(opt$out_prefix, "_genes-isoforms-reads_by_supercategory.svg"))



##---------------------------------------------------------
## STATS BREAKDOWN
##---------------------------------------------------------

novel_isof = c('Novel gene, antisense', 'Novel gene, divergent','Novel gene, intergenic','Novel gene, sense-overlap','Novel, monoexonic','Novel, in catalog','Novel, not in catalog')
novel_gene = c('Novel gene, antisense', 'Novel gene, divergent','Novel gene, intergenic','Novel gene, sense-overlap')

# Get some general stats
total_gene_count     = length( unique(d$gene_name) )
novel_gene_count     = length( unique(d$gene_name[d$supercategory %in% novel_gene]) )
novel_gene_pct       = round( (novel_gene_count / total_gene_count)*100, 0)

total_isof_count     = sum(d$isocount_all)
novel_isof_count     = sum(d$isocount_all[d$supercategory %in% novel_isof])
novel_isof_pct       = round( (novel_isof_count / total_isof_count)*100, 0)

total_read_count     = sum(d$readcount_all)
novel_isof_readcount = sum(d$readcount_all[d$supercategory %in% novel_isof])
novel_isof_readpct   = round( (novel_isof_readcount / total_read_count)*100, 0)

stats_out = data.frame(Category = c("total_gene_count", "novel_gene_count", "novel_gene_pct",
                                    "total_isof_count", "novel_isof_count", "novel_isof_pct",
                                    "total_read_count", "novel_isof_readcount", "novel_isof_readpct"),
                       Value    = c(total_gene_count, novel_gene_count, novel_gene_pct,
                                    total_isof_count, novel_isof_count, novel_isof_pct,
                                    total_read_count, novel_isof_readcount, novel_isof_readpct))
stats_out


# Write stats to file
write.table(stats_out, file=paste0(opt$out_prefix, "_gene-isoform-read_stats.txt"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


##-------------------------------------------------
## ISOFORM AND READ BREAKDOWNS ACROSS GENE BIOTYPES
##-------------------------------------------------

pd_iso  = d[, .(genecount = uniqueN(gene_name), isocount = sum(isocount_all), readcount = sum(readcount_all)), by = unified_type]

# For horizontal bars, reverse the factors
pd_iso$unified_type = factor(pd_iso$unified_type, levels = rev(levels(pd_iso$unified_type)))

# For horizontal bars, reverse the color palette (select the default ggplot colors from scales::hue_pal
default_colors  = scales::hue_pal()( length( levels(pd_iso$unified_type) ) )
reversed_colors = rev(default_colors)

# Calculate percentages to show next to the bars
short_label <- scales::label_number(scale_cut = cut_short_scale(), accuracy = 0.1)
pd_iso <- pd_iso %>%
    mutate(
        genepct   = genecount / sum(genecount),
        genelabel = paste0(short_label(genecount), " (", scales::percent(genepct, accuracy = 1), ")"),
        isopct    = isocount / sum(isocount),
        isolabel  = paste0(short_label(isocount),  " (", scales::percent(isopct, accuracy = 1),  ")"),
        readpct   = readcount / sum(readcount),
        readlabel = paste0(short_label(readcount), " (", scales::percent(readpct, accuracy = 1),  ")")
  )

# General plot elements
bar_width  = 0.8
text_size  = 24
label_size = 5

# Plot for gene count
plot_gcount = ggplot(pd_iso, aes(y = unified_type, x = genecount, fill = unified_type)) +
    geom_col( width = bar_width ) +
    geom_text( aes(label = genelabel), hjust = -0.1, size=label_size) +
    theme_minimal() +
    theme( text        = element_text(size=text_size),
           axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                     margin = margin(t = 10, b = 10, unit = "pt")),
           axis.text.y = element_text(hjust = 0),
           legend.position = "none") +
    labs(title = "Genes", x = "# Genes", y = "") +
    scale_fill_manual(values = reversed_colors) +
    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale(), accuracy = 1),
                       expand = expansion(mult = c(0, 0.5)))

# Plot for isoform count
plot_icount = ggplot(pd_iso, aes(y = unified_type, x = isocount, fill = unified_type)) +
    geom_col( width=bar_width ) +
    geom_text( aes(label = isolabel), hjust = -0.1, size=label_size) +
    theme_minimal() +
    theme( text            = element_text(size=text_size),
           legend.position = "none",
           axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                     margin = margin(t = 10, b = 10, unit = "pt")),
           axis.text.y     = element_blank()) +
    labs(title = "Isoforms", x = "# Isoforms", y = "") +
    scale_fill_manual(values = reversed_colors) +
    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale(), accuracy = 1),
                       expand = expansion(mult = c(0, 0.55)))

# Plot for read count
plot_rcount = ggplot(pd_iso, aes(y = unified_type, x = readcount, fill = unified_type)) +
    geom_col( width = bar_width ) +
    geom_text( aes(label = readlabel), hjust = -0.1, size=label_size) +
    theme_minimal() +
    theme( text            = element_text(size=text_size),
           legend.position = "none",
           axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                     margin = margin(t = 10, b = 10, unit = "pt")),
           axis.text.y     = element_blank()) +
    labs(title = "Reads", x = "# Reads", y = "", fill="Category") +
    scale_fill_manual(values = reversed_colors) +
    scale_x_continuous(labels = scales::label_number(scale_cut = cut_short_scale(), accuracy = 1),
                       expand = expansion(mult = c(0, 0.7)))

options(repr.plot.width=20, repr.plot.height=6)
iso_biotype = (plot_gcount | plot_icount | plot_rcount)
ggsave(plot=iso_biotype, file=paste0(opt$out_prefix, "_genes-isoforms-reads_by_biotype.pdf"), width=17, height=7)
ggsave(plot=iso_biotype, file=paste0(opt$out_prefix, "_genes-isoforms-reads_by_biotype.svg"), width=17, height=7)
clean_svg_file(paste0(opt$out_prefix, "_genes-isoforms-reads_by_biotype.svg"))



##-----------------------------------------------------
## PLOT ISOFORM DISTRIBUTIONS BY RANKED GENE EXPRESSION
##-----------------------------------------------------

pd_list = list()

for ( tg in tglist ){
    ## Make subset vector for this track group
    pd_sub = d[[paste0("isocount_", tg)]] > 0
    
    ## Select gene set
    pd_gene = d[ pd_sub, .(isocount = sum(isocount_all), readcount = sum(readcount_all)), by = gene_name]
    pd_gene = pd_gene[ order(pd_gene$readcount, decreasing=FALSE) ,]
    pd_gene$isocount  = -pd_gene$isocount
    pd_gene$readcount = log10(pd_gene$readcount)
    pd_gene$order = 1:nrow(pd_gene)
    head(pd_gene)

    # Create mirrored horizontal bar plot
    pd_list[[ paste0("iso_", tg) ]] = ggplot(pd_gene, aes(y = isocount, x=order)) +
        geom_col(fill="#00203FFF", width=1, color="#00203FFF" ) +
        coord_flip() +
        theme_minimal() +
        theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, 
                                         margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_blank(),
               plot.margin = margin(t = 5, r = 0, b = 5, l = 5, unit = "pt")) +
        labs(title = tg, x = "Expression Rank", y = "# Isoforms", fill="Category") +
        scale_y_continuous(labels = function(x) abs(x))

    pd_list[[ paste0("exp_", tg) ]] = ggplot(pd_gene, aes(y = readcount, x=order)) +
        geom_col(fill="#ADEFD1FF", width=1, color="#ADEFD1FF" ) +
        coord_flip() +
        theme_minimal() +
        theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x     = element_text(angle = 0, hjust = 1, vjust = 1, 
                                              margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_blank(),
               plot.margin     = margin(t = 5, r = 5, b = 5, l = 0, unit = "pt"),
               axis.title.y    = element_blank()
             ) +
        labs(title = "", x = "", y = "log10(Expression)", fill="Category")
}

# Output plots
options(repr.plot.width=18, repr.plot.height=12)
iso_ranking = (pd_list$iso_all | pd_list$exp_all)

ggsave(plot=iso_ranking, file=paste0(opt$out_prefix, "_isoform-expression-ranking.pdf"), width=8, height=8)
ggsave(plot=iso_ranking, file=paste0(opt$out_prefix, "_isoform-expression-ranking.svg"), width=8, height=8)
clean_svg_file(paste0(opt$out_prefix, "_isoform-expression-ranking.svg"))



##--------------------------------------------------------------
## ISOFORMS PER GENE, ISOFORM LENGTH, EXON, CODING, DISTRIBUTION
##--------------------------------------------------------------

## Gather datasets
pd_iso = d[ , .(isocount = sum(isocount_all)), by = c("gene_name", "supercategory")]
pd_len = d[ , c("supercategory", "length") ]
pd_exn = d[ , c("supercategory", "exon_count") ]
pd_nmd = d[ d$cds_type!="Non-coding" , .(total = .N, count_nmd = sum(cds_type == "NMD")), by = c("gene_name", "supercategory")][
             , pct_nmd := (count_nmd / total) * 100][
             , .(gene_name, supercategory, pct_nmd)]
pd_cds = d[ , .(total = .N, count_nmd = sum(cds_type != "Non-coding")), by = c("gene_name", "supercategory")][
             , pct_cds := (count_nmd / total) * 100][
             , .(gene_name, supercategory, pct_cds)]

# For horizontal bars, reverse the factors
pd_iso$supercategory = factor(pd_iso$supercategory, levels = rev(levels(d$supercategory)))
pd_len$supercategory = factor(pd_len$supercategory, levels = rev(levels(d$supercategory)))
pd_exn$supercategory = factor(pd_exn$supercategory, levels = rev(levels(d$supercategory)))
pd_nmd$supercategory = factor(pd_nmd$supercategory, levels = rev(levels(d$supercategory)))
pd_cds$supercategory = factor(pd_cds$supercategory, levels = rev(levels(d$supercategory)))

# Rescale
pd_len = pd_len[ pd_len$length < 10000, ]
pd_len$length = pd_len$length/1000

# For horizontal bars, reverse the color palette (select the default ggplot colors from scales::hue_pal
default_colors  = scales::hue_pal()( length( levels(pd_iso$supercategory) ) )
reversed_colors = rev(default_colors)

## Gather violin plots
boxplot_iso = ggplot(pd_iso, aes(x = supercategory, y = isocount, fill=supercategory)) +
  geom_boxplot(width=0.8, outlier.shape = NA ) +
  coord_flip() +
  labs(title = "", x = "", y = "# Isoforms") +
  theme_minimal() +
  theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x     = element_text(angle = 0, hjust = 1, vjust = 1, 
                                              margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_text(hjust = 0),
               plot.margin     = margin(r = 20, unit = "pt"),
             ) +
  ylim(0,15)+
  scale_fill_manual(values = reversed_colors)

boxplot_len = ggplot(pd_len, aes(x = supercategory, y = length, fill=supercategory)) +
  geom_boxplot(width=0.8, outlier.shape = NA ) +
  coord_flip() +
  labs(title = "", x = "", y = "Length (kb)") +
  theme_minimal() +
  theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x     = element_text(angle = 0, hjust = 1, vjust = 1, 
                                              margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_blank(),
               axis.title.y    = element_blank(), 
               plot.margin     = margin(r = 20, unit = "pt"),
             ) +
  scale_fill_manual(values = reversed_colors)

boxplot_exn = ggplot(pd_exn, aes(x = supercategory, y = exon_count, fill=supercategory)) +
  geom_boxplot(width=0.8, outlier.shape = NA ) +
  coord_flip() +
  labs(title = "", x = "", y = "# Exons") +
  theme_minimal() +
  theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x     = element_text(angle = 0, hjust = 1, vjust = 1, 
                                              margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_blank(),
               axis.title.y    = element_blank(), 
               plot.margin     = margin(r = 20, unit = "pt"),
             ) +
  ylim(0,30)+
  scale_fill_manual(values = reversed_colors)

boxplot_cds = ggplot(pd_cds, aes(x = supercategory, y = pct_cds, fill=supercategory)) +
  geom_boxplot(width=0.8, outlier.shape = NA ) +
  geom_jitter(width = 0.1, alpha = 0.5, size=0.3 ) +
  coord_flip() +
  labs(title = "", x = "", y = "% Coding") +
  theme_minimal() +
  theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x     = element_text(angle = 0, hjust = 1, vjust = 1, 
                                              margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_blank(),
               axis.title.y    = element_blank(), 
               plot.margin     = margin(r = 20, unit = "pt"),
             ) +
  scale_fill_manual(values = reversed_colors)

boxplot_nmd = ggplot(pd_nmd, aes(x = supercategory, y = pct_nmd, fill=supercategory)) +
  geom_boxplot(width=0.8, outlier.shape = NA ) +
  geom_jitter(width = 0.1, alpha = 0.5, size=0.3 ) +
  coord_flip() +
  labs(title = "", x = "", y = "% NMD") +
  theme_minimal() +
  theme( text            = element_text(size=text_size),
               legend.position = "none",
               axis.text.x     = element_text(angle = 0, hjust = 1, vjust = 1, 
                                              margin = margin(t = 10, b = 10, unit = "pt")),
               axis.text.y     = element_blank(),
               axis.title.y    = element_blank(), 
               plot.margin     = margin(r = 20, unit = "pt"),
             ) +
  scale_fill_manual(values = reversed_colors)

options(repr.plot.width=20, repr.plot.height=8)
iso_stats = boxplot_iso | boxplot_len | boxplot_exn | boxplot_cds | boxplot_nmd

ggsave(plot=iso_stats, file=paste0(opt$out_prefix, "_isoform-stats-breakdown.pdf"), width=20, height=8)
ggsave(plot=iso_stats, file=paste0(opt$out_prefix, "_isoform-stats-breakdown.svg"), width=20, height=8)
clean_svg_file(paste0(opt$out_prefix, "_isoform-stats-breakdown.svg"))


# CLEANUP: Remove the empty Rplots.pdf generated by aplot/ggbreak alignment checks
if (file.exists("Rplots.pdf")) {
    file.remove("Rplots.pdf")
}
