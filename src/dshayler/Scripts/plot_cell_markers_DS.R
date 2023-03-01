#!/usr/binRscript

library(stchlk.heatmaply)
library(biomaRt)
library(shiny)
library(heatmaply)
library(shinyHeatmaply)

pdf(NULL)

# define mart for biomaRt gene lookup -------------------------------------
ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")


# load gene_transcript assocication table ---------------------------------
t_to_g_association <- readRDS("~/Homo_sapiens/grch38_tran/t_to_g_assocation.rds")


# define markers of intereste ---------------------------------------------
marker_list <- list(cone_markers     = c("OPN1LW", "RXRG", "CRX", "ARR3", "RORB", "GTF2IRD1", "GNAT2", "THRB", "MDM2", "MYCN"),
                    rod_markers      = c("NRL", "NXNL1"),
                    ama_hoz_markers  = c("PROX1", "PAX6", "NES"),
                    prog_mul_markers = c("NES", "RLBP1", "SOX2"))

# load datasets (census matrix) -------------------------------------------

dataset_paths <- c(
  dshayler_20170407_census = "~/single_cell_pipeline/output/FACS_20170407_dshayler_H_sapiens_output/Dominik_stringtie.tpm_census_matrix.csv",
  dshayler_20171031_census = "~/single_cell_pipeline/output/FACS_20171031_dshayler_H_sapiens_output/stringtie_transcripts.tpm_census_matrix.csv")

datasets <- lapply(dataset_paths, read.table, sep = "\t", header = TRUE)

# create a list of gene_count matrices on marker sets -----------------------------
gene_counts <- lapply(datasets, return_counts_par, marker_list)

#' save gene count matrices as .rda files for easy loading into heatmap
save_paths <- c(
  dshayler_20170407_cell_markers = "~/single_cell_pipeline/results/dshayler/marker_gene_count_mat_20170407_dshayler.rda",
  dshayler_20171031_cell_markers = "~/single_cell_pipeline/results/dshayler/marker_gene_count_mat_20171031_dshayler.rda")

purrr::map2(gene_counts, save_paths, function(x,y) {save(x, file = y)})




#' run heatmap visualization in shiny
runApp(system.file("shinyapp", package = "shinyHeatmaply"))

View(gene_counts)
