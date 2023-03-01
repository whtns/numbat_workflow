## ----setup, echo = F-----------------------------------------------------
knitr::opts_chunk$set(echo = F, message = F, warning = F)


## ------------------------------------------------------------------------
library('tidyverse')
 library('fs')
 library('rprojroot')
library(seuratTools)
proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))


## ------------------------------------------------------------------------
txi_genes <- seuratTools::load_counts_from_stringtie(proj_dir, txOut = F)
txi_transcripts <- seuratTools::load_counts_from_stringtie(proj_dir, txOut = T)


## ------------------------------------------------------------------------
tpm_meta <- seuratTools::load_meta(proj_dir)


## ------------------------------------------------------------------------
feature_seus <- map(list(gene = txi_genes, transcript = txi_transcripts), seu_from_tibbles, tpm_meta)


## ------------------------------------------------------------------------
feature_seus_thresh <- map(feature_seus, filter_low_rc_cells)


## ------------------------------------------------------------------------
purrr::imap(feature_seus_thresh, save_seurat)

