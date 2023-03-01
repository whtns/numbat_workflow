#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(clustifyr)
library(Seurat)
library(patchwork)
library(numbat)
library(glue)

# work with numbat object ------------------------------

save_numbat_as_rds <- function(numbat_dir){
  # browser()
  nb = Numbat$new(out_dir = numbat_dir)

  numbat_rds_path <- paste0(numbat_dir, "_numbat.rds")

  saveRDS(nb, numbat_rds_path)

  return(nb)

}

safe_nasrds <- safely(save_numbat_as_rds, otherwise = NA_real_)

plot_phylo_w_celltypes <- function(seu_path, nb) {
  # browser()
  myseu <- readRDS(seu_path)

  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    # as.vector() %>%
    tibble::deframe() %>%
    set_names(str_replace(names(.), "\\.", "-")) %>%
    identity()

  mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  nb$plot_phylo_heatmap(
    clone_bar = TRUE,
    p_min = 0.9,
    pal_clone = mypal,
    annot = celltypes
  ) +
    labs(title = fs::path_file(seu_path))
}

safe_plot_phylo <- safely(plot_phylo_w_celltypes, otherwise = NA_real_)

# field ------------------------------
study = "field"

numbat_dirs <-
  dir_ls(glue("output/numbat/{study}_et_al/"), regexp = "/SRR[0-9]*$") %>%
  fs::path_abs() %>%
  set_names(fs::path_file(.))

seu_paths <- dir_ls(glue("output/seurat/{study}_et_al"), regexp = "SRR[0-9]*_infercnv_numbat_seu.rds")

nbs <- purrr::map(numbat_dirs, safe_nasrds)

nbs0 <- map(nbs, "result")

phylo_heatmaps <- purrr::map2(seu_paths, nbs0[c(2:5)], safe_plot_phylo)

pdf(glue("results/{study}/phylo_heatmaps.pdf"))
phylo_heatmaps
dev.off()

# collin ------------------------------
study = "collin"

numbat_dirs <-
  dir_ls(glue("output/numbat/{study}_et_al/"), regexp = "/SRR[0-9]*$") %>%
  fs::path_abs() %>%
  set_names(fs::path_file(.))

seu_paths <- dir_ls(glue("output/seurat/{study}_et_al"), regexp = "SRR[0-9]*_infercnv_numbat_seu.rds")

nbs <- purrr::map(numbat_dirs, safe_nasrds)

nbs0 <- map(nbs, "result")

phylo_heatmaps <- purrr::map2(seu_paths, nbs0[c(2:5)], safe_plot_phylo)

pdf(glue("results/{study}/phylo_heatmaps.pdf"))
phylo_heatmaps
dev.off()

# yang ------------------------------
study = "yang"

numbat_dirs <-
  dir_ls(glue("output/numbat/{study}_et_al/"), regexp = "/SRR[0-9]*$") %>%
  fs::path_abs() %>%
  set_names(fs::path_file(.))

seu_paths <- dir_ls(glue("output/seurat/{study}_et_al"), regexp = "SRR[0-9]*_infercnv_numbat_seu.rds")

# nbs <- purrr::map(numbat_dirs, safe_nasrds)
# nbs0 <- map(nbs, "result")

nbs0 <- dir_ls(glue("output/numbat/{study}_et_al"), regexp = "SRR[0-9]*_numbat.rds") %>%
  map(readRDS)

phylo_heatmaps <- purrr::map2(seu_paths, nbs0[c(2:5)], safe_plot_phylo)

pdf(glue("results/{study}/phylo_heatmaps.pdf"))
phylo_heatmaps
dev.off()

# wu ------------------------------
study = "wu"

numbat_dirs <-
  dir_ls(glue("output/numbat/{study}_et_al/"), regexp = "/SRR[0-9]*$") %>%
  fs::path_abs() %>%
  set_names(fs::path_file(.))

seu_paths <- dir_ls(glue("output/seurat/{study}_et_al"), regexp = "SRR[0-9]*.*_seu.rds") %>%
  str_subset(pattern = ".*dropped.*", negate = TRUE) %>%
  set_names(path_file(.))

nbs <- purrr::map(numbat_dirs, safe_nasrds)
nbs0 <- map(nbs, "result")

# nbs0 <- dir_ls(glue("output/numbat/{study}_et_al"), regexp = "SRR[0-9]*_numbat.rds") %>%
# map(readRDS) %>%
# identity()

phylo_heatmaps <- purrr::map2(seu_paths, nbs0, safe_plot_phylo)

pdf(glue("results/{study}/phylo_heatmaps.pdf"))
phylo_heatmaps
dev.off()

browseURL("results/wu/phylo_heatmaps.pdf")

#
# subset_numbat_object <- function(nb, drop_cells){
#
#   # rownames------------------------------
#   c("gexp_roll_wide")
#
#   # cell------------------------------
#   list_cats <- c("allele_post", "clone_post", "joint_post")
#
#   mynb <- nb
#
#   for (i in list_cats){
#     mynb[[i]] <- mynb[[i]][!mynb[[i]]$cell %in% drop_cells,]
#   }
#
#   mynb[["gtree"]] <-
#     mynb[["gtree"]] %>%
#     tidygraph::activate(nodes) %>%
#     filter(!name %in% drop_cells) %>%
#     identity()
#
#   return(mynb)
# }
#
# rod_cells <- str_replace(WhichCells(myseu, expression = type == "Rods"), "\\.", "-")
#
# cone_cells <- str_replace(WhichCells(myseu, expression = type == "Cones"), "\\.", "-")
#
# nb_subset <- nb %>%
#   subset_numbat_object(cone_cells)
#
# mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")
#
# # debug(numbat:::plot_phylo_heatmap)
# nb_subset$plot_phylo_heatmap(
#   clone_bar = TRUE,
#   p_min = 0.5,
#   pal_clone = mypal
# )





