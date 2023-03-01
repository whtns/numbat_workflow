#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
library(glue)

plot_phylo_w_celltypes <- function(nb, seu_path) {
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
    annot = celltypes,
    show_phylo = FALSE
  ) +
    labs(title = fs::path_file(seu_path))
}

safe_plot_phylo <- safely(plot_phylo_w_celltypes, otherwise = NA_real_)

# wu ------------------------------

study = "wu"

nb_paths <- dir_ls(glue("output/numbat/{study}_et_al/"), glob = "*_numbat.rds") %>%
  set_names(path_file(.)) %>%
  identity()

nbs <- nb_paths %>%
  map(readRDS)

seu_paths <- dir_ls(glue("output/seurat/{study}_et_al/")) %>%
  str_subset("dropped", negate = TRUE) %>%
  identity()

seu_paths <- seu_paths[c(2,4,5,8,9)]

phylo_heatmaps <- map2(nbs, seu_paths, safe_plot_phylo)

phylo_pdf <- glue("results/{study}/phylo_heatmaps.pdf")

pdf(phylo_pdf)
phylo_heatmaps
dev.off()

browseURL(phylo_pdf)


heatmap_file <- glue("results/{study}/phylo_heatmap_extra.pdf")

pdf(heatmap_file)
phylo_heatmap
dev.off()

browseURL(heatmap_file)


# yang ------------------------------

study = "yang"

nb_paths <- dir_ls(glue("output/numbat/{study}_et_al/"), glob = "*_numbat.rds") %>%
  set_names(path_file(.)) %>%
  identity()

nbs <- nb_paths %>%
  map(readRDS)

seu_paths <- dir_ls(glue("output/seurat/{study}_et_al/")) %>%
  str_subset("dropped", negate = TRUE) %>%
  identity()

seu_paths <- seu_paths[c(2,4,5,8,9)]

phylo_heatmaps <- map2(nbs, seu_paths, safe_plot_phylo)

phylo_pdf <- glue("results/{study}/phylo_heatmaps.pdf")

pdf(phylo_pdf)
phylo_heatmaps
dev.off()

browseURL(phylo_pdf)


heatmap_file <- glue("results/{study}/phylo_heatmap_extra.pdf")

pdf(heatmap_file)
phylo_heatmap
dev.off()

browseURL(heatmap_file)
