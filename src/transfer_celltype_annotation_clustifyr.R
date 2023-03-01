#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(clustifyr)
library(Seurat)

collin_ref <- readRDS("~/Homo_sapiens/numbat/collin_ref.rds")

sridhar_ref <- readRDS("~/Homo_sapiens/numbat/sridhar_ref.rds")

plae_ref <- readRDS("~/Homo_sapiens/numbat/plae_ref.rds")


drop_bad_cells <- function(seu_path, bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))) {
  browser()
  seu <- readRDS(seu_path)

  seu <- seu[,!seu$type %in% bad_cell_types]

  saveRDS(seu, str_replace(seu_path, "_seu.rds", "_dropped_cells_seu.rds"))

  retainedcells_dimplot <- DimPlot(
    seu,
    group.by = "type"
  ) +
    labs(title = fs::path_file(seu_path))

  return(seu_path)

}

plot_celltype_predictions <- function(seu_path, bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))) {
  browser()
  seu <- readRDS(seu_path)

  seu_mat <- GetAssayData(seu, slot = "data")

  res <- clustify(
    input = seu_mat,
    metadata = seu$gene_snn_res.0.2,
    ref_mat = plae_ref,
    query_genes = VariableFeatures(seu)
  )

  cor_to_call(res)

  res2 <- cor_to_call(
    cor_mat = res,                  # matrix correlation coefficients
    cluster_col = "gene_snn_res.0.2" # name of column in meta.data containing cell clusters
  )

  # Insert into original metadata as "type" column
  seu@meta.data <- call_to_metadata(
    res = res2,                     # data.frame of called cell type for each cluster
    metadata = seu@meta.data,           # original meta.data table containing cell clusters
    cluster_col = "gene_snn_res.0.2" # name of column in meta.data containing cell clusters
  )

  saveRDS(seu, seu_path)

  allcells_dimplot <- DimPlot(
    seu,
    group.by = "type"
  ) +
    labs(title = fs::path_file(seu_path))

  seu <- seu[,!seu$type %in% bad_cell_types]

  saveRDS(seu, str_replace(seu_path, "_seu.rds", "_dropped_cells_seu.rds"))

  retainedcells_dimplot <- DimPlot(
    seu,
    group.by = "type"
  ) +
    labs(title = fs::path_file(seu_path))

  return(list("allcells" = allcells_dimplot, "retainedcells" = retainedcells_dimplot))

}

transfer_celltype_predictions <- function(seu_path, bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))) {
  # browser()
  seu <- readRDS(seu_path)

  numbat_seu_path = str_replace(seu_path, "_seu.rds", "_cnv_seu.rds")

  numbat_seu <- readRDS(numbat_seu_path)

  numbat_seu <- AddMetaData(numbat_seu, seu@meta.data[c("type", "r")])

  saveRDS(numbat_seu, numbat_seu_path)

  return(numbat_seu_path)
}

# field ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/field_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_seu.rds")

celltype_predictions <-
  seu_paths %>%
  map(plot_celltype_predictions)

pdf_path <- "~/single_cell_projects/resources/field_et_al_proj/results/field_celltype_predictions.pdf"

pdf(pdf_path)
celltype_predictions
dev.off()

browseURL(pdf_path)

# collin ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/collin_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_seu.rds")

celltype_predictions <-
  seu_paths %>%
  map(plot_celltype_predictions)

pdf_path <- "~/single_cell_projects/resources/collin_et_al_proj/results/collin_celltype_predictions.pdf"

pdf(pdf_path)
celltype_predictions
dev.off()

browseURL(pdf_path)

# yang ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/yang_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_seu.rds")

celltype_predictions <-
  seu_paths %>%
  map(plot_celltype_predictions)

pdf_path <- "~/single_cell_projects/resources/yang_et_al_proj/results/yang_celltype_predictions.pdf"

pdf(pdf_path)
celltype_predictions
dev.off()

browseURL(pdf_path)

# wu ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/wu_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_seu.rds")

celltype_predictions <-
  seu_paths %>%
  map(plot_celltype_predictions)

pdf_path <- "~/single_cell_projects/resources/wu_et_al_proj/results/wu_celltype_predictions.pdf"

pdf(pdf_path)
celltype_predictions
dev.off()

browseURL(pdf_path)
