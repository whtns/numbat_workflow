#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(clustifyr)
library(Seurat)
library(patchwork)
library(numbat)
library(glue)

collin_ref <- readRDS("~/Homo_sapiens/numbat/collin_ref.rds")

sridhar_ref <- readRDS("~/Homo_sapiens/numbat/sridhar_ref.rds")

plae_ref <- readRDS("~/Homo_sapiens/numbat/plae_ref.rds")

calc_celltype_proportion <- function(seu_path) {
  # browser()
  seu <- readRDS(seu_path)

  seu_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell")

  return(seu_meta)

}

#define custom color scale
present_cell_types <- c("Rods", "Muller Glia", "Microglia", "Cones", "RPE", "Red Blood Cells",
                        "Neurogenic Cells", "Endothelial", "Schwann", "Fibroblasts",
                        "Pericytes", "RPCs", "Astrocytes", "Bipolar Cells", "Horizontal Cells",
                        "Rod Bipolar Cells", "Late RPCs")

myColors <- scales::hue_pal()(length(present_cell_types))
names(myColors) <- present_cell_types
custom_colors <- scale_fill_manual(name = "type", values = myColors)

make_celltype_proportion_plot <- function(celltype_proportion_df, custom_colors){
  proportion_plot <- celltype_proportion_df %>%
    dplyr::mutate(sample_id = str_remove(sample_id, "_seu.rds")) %>%
    ggplot(aes(sample_id, fill = type)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    custom_colors +
    coord_flip()
}


# field ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/field_et_al_proj/output/seurat/", regexp = ".*_seu.rds") %>%
  str_subset(pattern = ".*dropped.*", negate = TRUE)

field_celltype_proportions <-
  seu_paths %>%
  set_names(path_file(.)) %>%
  map_dfr(calc_celltype_proportion, .id = "sample_id")

# plot_celltype_proportions

field_proportions <-
  make_celltype_proportion_plot(field_celltype_proportions, custom_colors) +
  labs(title = "Field et al.")

# collin ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/collin_et_al_proj/output/seurat/", regexp = "*_seu.rds") %>%
  str_subset(pattern = ".*dropped.*", negate = TRUE)

collin_celltype_proportions <-
  seu_paths %>%
  set_names(path_file(.)) %>%
  map_dfr(calc_celltype_proportion, .id = "sample_id")

# plot_celltype_proportions

collin_proportions <-
  make_celltype_proportion_plot(collin_celltype_proportions, custom_colors) +
  labs(title = "Collin et al.")

# yang ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/yang_et_al_proj/output/seurat/", regexp = "*_seu.rds") %>%
  str_subset(pattern = ".*dropped.*", negate = TRUE)

yang_celltype_proportions <-
  seu_paths %>%
  set_names(path_file(.)) %>%
  map_dfr(calc_celltype_proportion, .id = "sample_id")

# plot_celltype_proportions

yang_proportions <-
  make_celltype_proportion_plot(yang_celltype_proportions, custom_colors) +
  labs(title = "Yang et al.")

# wu ------------------------------
seu_paths <- dir_ls("~/single_cell_projects/resources/wu_et_al_proj/output/seurat/", regexp = "*_seu.rds") %>%
  str_subset(pattern = ".*dropped.*", negate = TRUE)

wu_celltype_proportions <-
  seu_paths %>%
  set_names(path_file(.)) %>%
  map(calc_celltype_proportion)

wu_celltype_proportions$SRR13884246_infercnv_numbat_seu.rds$GT_opt <- as.character(wu_celltype_proportions$SRR13884246_infercnv_numbat_seu.rds$GT_opt)

wu_celltype_proportions <-
  wu_celltype_proportions %>%
  # map(mutate, GT_opt = as.character(GT_opt)) %>%
  dplyr::bind_rows(.id = "sample_id") %>%
  identity()

# plot_celltype_proportions

wu_proportions <-
  make_celltype_proportion_plot(wu_celltype_proportions, custom_colors) +
  labs(title = "Wu et al.")

celltypes_by_study_plot <- collin_proportions + field_proportions + wu_proportions + yang_proportions + plot_layout(guides = "collect")

metadata_by_study_table <-
  list("collin" = collin_celltype_proportions,
       "wu" = wu_celltype_proportions,
       "field" = field_celltype_proportions,
       "yang" = yang_celltype_proportions) %>%
  dplyr::bind_rows(.id = "study") %>%
  dplyr::mutate(type = dplyr::coalesce(type, type.clustify)) %>%
  dplyr::mutate(r = dplyr::coalesce(r, r.clustify)) %>%
  dplyr::select(-type.clustify, -r.clustify)

write_csv(metadata_by_study_table, "results/metadata_by_study.csv")

celltype_proportions_by_sample <-
  metadata_by_study_table %>%
  dplyr::group_by(sample_id, type) %>%
  dplyr::count(type) %>%
  tidyr::pivot_wider(names_from = "type", values_from = "n") %>%
  mutate(across(where(is.integer), replace_na, 0)) %>%
  tidyr::pivot_longer(-sample_id, names_to = "type", values_to = "n") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(percent_by_sample = n/sum(n)) %>%
  group_by(type) %>%
  summarize(mean_percent = mean(percent_by_sample), mean_n = mean(n)) %>%
  dplyr::arrange(desc(mean_percent)) %>%
  identity()

write_csv(celltype_proportions_by_sample, "results/celltypes_by_sample.csv")

test0 <- janitor::tabyl(metadata_by_study_table, type, clone_opt, sample_id)




