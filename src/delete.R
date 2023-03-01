#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(SeuratDisk)


external_sc_rb_projects <- c(collin = "~/single_cell_projects/resources/collin_et_al_proj/output/seurat/",
  wu = "~/single_cell_projects/resources/wu_et_al_proj/output/seurat/",
  field = "~/single_cell_projects/resources/field_et_al_proj/output/seurat/",
  yang = "~/single_cell_projects/resources/yang_et_al_proj/output/seurat/")

seu_paths <- map(external_sc_rb_projects, ~dir_ls(.x, regexp = "SRR[0-9]*_seu.rds", recurse = TRUE)) %>%
  unlist()

convert_to_seuratdisk <- function(seu_path){
  seu <- readRDS(seu_path)

  seuratdisk_path <- str_replace(seu_path, ".rds", ".h5Seurat")

  SaveH5Seurat(seu, filename = seuratdisk_path, overwrite = TRUE)
}

convert_to_seuratdisk(seu_paths[[1]])


