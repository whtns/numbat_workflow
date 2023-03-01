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
  # browser()
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

# field ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/field_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_infercnv_numbat_seu.rds")

dropped_seus <-
  seu_paths %>%
  map(drop_bad_cells)

# collin ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/collin_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_seu.rds")

dropped_seus <-
  seu_paths %>%
  map(drop_bad_cells)


# yang ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/yang_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_infercnv_numbat_seu.rds")

dropped_seus <-
  seu_paths %>%
  map(drop_bad_cells)

# wu ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/wu_et_al_proj/output/seurat/", regexp = "SRR[0-9]*_infercnv_numbat_seu.rds")

dropped_seus <-
  seu_paths %>%
  map(drop_bad_cells)
