#!/usr/bin/Rscript

#load required libraries and functions
#=====================================
library(optparse)

find_cells <- function(df, day_cat, treatment){
  df <- dplyr::filter(df, day %in% day_cat) %>% 
    dplyr::filter(treatment_group %in% treatment) %>% 
    dplyr::filter(branch == "A")
  cells = paste(df$Sample_ID, collapse = " ")
  return(cells)
}

# SHL 20170407 ------------------------------------------------------------
# default_cell_info = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/cell_info.csv"
# default_feature = "transcript"
# default_groupfile = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_deseq2_comparison_file.csv"
# default_infile = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/transcript_count_matrix.csv"
# # default_infile = "/home/skevin/tmp/FACS_20170407_sunlee_census_matrix_mod.csv" # census
# # default_infile = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/gene_count_matrix.csv"
# default_out = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/"
# # raw transcript counts
# input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_differential_expression_input.rda"
# # census transcript counts
# # input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_census_differential_expression_input.rda"
# # raw gene counts
# # input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_genes_differential_expression_input.rda"

# DS 2seq ------------------------------------------------------------
# default_cell_info = "/home/skevin/tmp/2_seq_meta_111517.csv"
# # default_cell_info = "~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/2_seq_dshayler/2_seq_meta_111517.csv"
# default_feature = "transcript"
# default_groupfile = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_diffex_comparison_file.csv"
# default_infile = "/home/skevin/tmp/FACS_20170407_20171031_dshayler_transcripts_raw_counts.csv"
# # default_infile = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_transcripts_raw_counts.csv"
# default_out = "~/single_cell_pipeline/output/merged_analyses/"
# # # raw transcript counts
input_rda <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_differential_expression_input.rda"

# diffex_input <- mget(ls(pattern = "default"))
# save(diffex_input, file = input_rda)

if (file.exists(input_rda)){
  load( input_rda)
  list2env(diffex_input, globalenv())
} else {
  default_expr_mat = default_annotation =  default_cell_info =  default_plot_settings = default_out = NA
  
}


# default input files -----------------------------------------------------

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-i", "--indir"), type="character", default=default_out,
              help="heatmap input directory (same as pior deseq2 output) [default= %default]", metavar="character"),
  make_option(c("-g", "--groupfile"), type="character", default=default_groupfile,
              help="file that defines comparison groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=default_out,
              help="output directory name [default= %default]", metavar="character"), 
  make_option(c("-f", "--feature"), type="character", default=NA, 
              help="the feature (gene or transcript) under analysis", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  missing_args <- paste(names(which(sapply(opt, is.na))), collapse = " ")
  stop(paste0("Please provide all necessary arguments. Missing: ", missing_args), call.=FALSE)
}


# load required libraries -------------------------------------------------
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
suppressMessages(library(gtools))
suppressMessages(library(data.table))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(shiny))
suppressMessages(library(heatmaply))
suppressMessages(library(shinyHeatmaply))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
edb <- EnsDb.Hsapiens.v86

lookup_genes <- function(txname){
  
  txs <- transcripts(edb, filter = TxIdFilter(txname),
                     columns = c("symbol"))
  return(txs$symbol)
  
}

lookup_transcripts <- function(genename){
  
  txs <- transcripts(edb, filter = GenenameFilter(genename),
                     columns = c("tx_id"))
  return(txs$tx_id)
  
}

user_groups <- read.table(textConnection(gsub(" ", "\t", readLines(opt$groupfile))), sep="\t", header = TRUE, stringsAsFactors = FALSE, fill=TRUE)

# load deseq2 deseqdataset from LRT test

ddsRDS <- paste0(opt$indir, user_groups$analysis_id, "_dds.rds")
ddsl <- lapply(ddsRDS, readRDS)

resl <- purrr::map(ddsl, cataract::get_deseq2_res, feature = opt$feature)

# load deseq2 results of LRT test
resRDS <- paste0(opt$indir, user_groups$analysis_id, "_res.rds")
purrr::map2(resl, resRDS, function(x,y) saveRDS(x, file = y))

resl <- lapply(resRDS, readRDS)

# load vst transform of the data set
vsdRDS <- paste0(opt$indir, user_groups$analysis_id, "_vsd.rds")

vsdl <- lapply(ddsl, vst, blind = F)
purrr::map2(vsdl, vsdRDS, function(x,y) saveRDS(x, file=y))

vsdl <- lapply(vsdRDS, readRDS)

# plot cell by gene heatmap -----------------------------------------------
find_top_var_genes <- function(vsd_obj, top_n){
  vsd <- head( order( rowVars( assay(vsd_obj) ), decreasing=TRUE ), top_n )
  vsd <- assay(vsd_obj)[vsd,]
  syms <- lookup_genes(rownames(vsd))
  rownames(vsd) <- paste0(syms, "\n", rownames(vsd))
  return(vsd)
} 

find_top_sig_genes <- function(vsd_obj, res, top_n){
  vsp <- rownames(head(res[order(res$padj),], top_n))
  vsp <- assay(vsd_obj)[vsp,]
  syms <- lookup_genes(rownames(vsp))
  rownames(vsp) <- paste0(syms, "\n", rownames(vsp))
  return(vsp)
} 

find_user_supplied_features <- function(vsd_obj, in_features, top_n){
  in_features <- read.csv(in_features, header = FALSE)
  vsp <- in_features[[1]]
  vsp <- assay(vsd_obj)[vsp,]
  syms <- lookup_genes(rownames(vsp))
  rownames(vsp) <- paste0(syms, "\n", rownames(vsp))
  return(vsp)
} 


topVarGenes <- lapply(vsdl, find_top_var_genes, 500)
names(topVarGenes) <- gsub("_res.rds", "_top_var", basename(resRDS))

topSigGenes <- purrr::map2(vsdl, resl, find_top_sig_genes, 500)
names(topVarGenes) <- gsub("_res.rds", "_top_sig", basename(resRDS))

user_genes <- "~/single_cell_tools/test_pcs_4_6_7.csv"
topUserGenes <- purrr::map(vsdl, find_user_supplied_features, user_genes, 500)
names(topUserGenes) <- gsub("_res.rds", "_top_user", basename(resRDS))

# heatmap2 method ---------------------------------------------------------
# hm <- heatmap.2( topVarGenes[[1]], scale="row",
#                  trace="none", dendrogram="column",
#                  col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

plot_heatmap <- function(filt_vsd, heatmap_path, file_tag){
  heatmap_path <- paste0(heatmap_path, file_tag, "_heatmap.html")
  heatmaply(filt_vsd, scale = "row", seriate = "none", dendrogram = "both", hclust_method = "ward", file = heatmap_path, colors = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255) )
  
}

#test_filtv <- assay(vsdl[[1]])[topVarGenes[[1]],]

#test_hm_path <- heatmap_paths[[1]]

heatmap_paths <- paste0(opt$out, user_groups$analysis_id)

purrr::map2(topSigGenes, heatmap_paths, plot_heatmap, "_top_sig_genes")
purrr::map2(topVarGenes, heatmap_paths, plot_heatmap, "_top_var_genes")
purrr::map2(topUserGenes, heatmap_paths, plot_heatmap, "_top_user_genes")
