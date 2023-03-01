#!/usr/bin/env Rscript

## ----setup, echo = F-----------------------------------------------------
knitr::opts_chunk$set(echo = F, message = F, warning = F)


## ------------------------------------------------------------------------
library(tidyverse)
library(velocyto.R)
library(Seurat)
library(rprojroot)
library(fs)
library(gsubfn)
proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

## ------------------------------------------------------------------------

set_loom_names <- function(loom_mats, batch){
  # browser()
  oldnames <- colnames(loom_mats[[1]])
  
  toreplace<-list("MyTissue:" = paste0(batch, "_S"), "_removed_duplicates.bam" = "")
  newnames <- gsubfn(paste(names(toreplace),collapse="|"),toreplace,oldnames)
  
  name_loom_mat <- function(loom_mat, newnames){
    colnames(loom_mat) <- newnames
    return(loom_mat)
  }
  
  named_loom_mats <- purrr::map(loom_mats, name_loom_mat, newnames)	
  return(named_loom_mats)
}

save_anno <- function(seu, path){
  anno_tbl <- as_tibble(seu[[]])
  write_csv(anno_tbl, path)
}



## ------------------------------------------------------------------------
### ASSUMES Seurat object has already been computed -- in this case the object tseu is the result of seurat

### load seurat object
features = c("gene", "transcript")
seu_paths <- fs::path(proj_dir, "output", "sce") %>% 
  dir_ls() %>% 
  path_filter("*_seu*.rds") %>%
  identity()

# seus <- map(seu_paths, readRDS)

# ------------------------------------------------------------------------
# IMPORTANT!!! create a loom file using the velocyto command line tools first
#line 114 cannot run in chunk; use terminal

loom_path <- fs::path(proj_dir, "output", "velocyto", "4_seq_dshayler_fetal_HS.loom") %>% 
  identity()

ldat <- read.loom.matrices(loom_path)

set_velocity_params <- function(seu_paths, filter_type, feature_type, embedding) {
  print(paste0("running ", paste0("fetal_", feature_type, filter_type)))
  # browser()  
  seu_path <- seu_paths %>% 
    path_filter(paste0("*", feature_type, "_seu", filter_type, ".rds"))
  
  seu <- readRDS(seu_path)
  
  velocity_params <- list(
    seu = seu,
    seu_path = seu_path,
    ldat = ldat,
    tt_flag = "fetal",
    feature_flag = feature_type,
    filter_flag = filter_type,
    emb_flag = embedding
  )
  
  return(velocity_params)
}

# gene unfiltered
velocity_params <- set_velocity_params(seu_paths, filter_type = "", feature_type = "gene", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# gene nonPRs
velocity_params <- set_velocity_params(seu_paths, filter_type = "_remove_nonPRs", feature_type = "gene", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# gene lowrc
velocity_params <- set_velocity_params(seu_paths, filter_type = "_remove_lowrc", feature_type = "gene", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# gene nonPRs and lowrc
velocity_params <- set_velocity_params(seu_paths, filter_type = "_remove_nonPRs_and_lowrc", feature_type = "gene", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# transcript unfiltered
velocity_params <- set_velocity_params(seu_paths, filter_type = "", feature_type = "transcript", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# transcript nonPRs
velocity_params <- set_velocity_params(seu_paths, filter_type = "_remove_nonPRs", feature_type = "transcript", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# transcript lowrc
velocity_params <- set_velocity_params(seu_paths, filter_type = "_remove_lowrc", feature_type = "transcript", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

# transcript nonPRs and lowrc
velocity_params <- set_velocity_params(seu_paths, filter_type = "_remove_nonPRs_and_lowrc", feature_type = "transcript", embedding = "umap")

source(fs::path(proj_dir, "src", "rna_velocity", "run_velocyto.R"))

