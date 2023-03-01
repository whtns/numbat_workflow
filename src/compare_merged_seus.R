#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(Seurat)
library(SeuratDisk)


Convert("~/single_cell_projects/resources/yang_et_al_proj/output/scanpy/merged/merged_study_processed.h5ad", dest = "h5seurat", overwrite = TRUE)


merged_seu <- LoadH5Seurat("~/single_cell_projects/resources/yang_et_al_proj/output/scanpy/merged/merged_study_processed.h5seurat")

