#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)
library(glue)


retrieve_cell_stats <- function(seu_path){

  seu <- readRDS(seu_path)

  stats = seu@meta.data[c("nCount_gene", "nFeature_gene", "percent.mt")] %>%
    tibble::rownames_to_column("cell")

  return(stats)
}

# collin ------------------------------

study = "collin"

collin_cell_stats <- dir_ls(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat"), regexp = "SRR[0-9]*_seu.rds") %>%
  set_names(path_file(.)) %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# field ------------------------------

study = "field"

field_cell_stats <- dir_ls(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat"), regexp = "SRR[0-9]*_seu.rds") %>%
  set_names(path_file(.)) %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# wu ------------------------------

study = "wu"

wu_cell_stats <- dir_ls(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat"), regexp = "SRR[0-9]*_seu.rds") %>%
  set_names(path_file(.)) %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# yang ------------------------------

study = "yang"

yang_cell_stats <- dir_ls(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat"), regexp = "SRR[0-9]*_seu.rds") %>%
  set_names(path_file(.)) %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

study_cell_stats <- dplyr::bind_rows(list("collin" = collin_cell_stats, "field" = field_cell_stats, "wu" = wu_cell_stats, "yang" = yang_cell_stats), .id = "study")

umis_per_cell <- ggplot(study_cell_stats, aes(x = nCount_gene, y = sample_id, fill = study, group = sample_id)) +
  geom_density_ridges() +
  scale_x_log10() +
  theme(
    axis.title.y=element_blank(),  #remove y axis labels,
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()  #remove y axis ticks
  ) +
  labs(title = "UMIs per cell")

genes_per_cell <- ggplot(study_cell_stats, aes(x = nFeature_gene, y = sample_id, fill = study, group = sample_id)) +
  geom_density_ridges() +
  scale_x_log10() +
  theme(
    axis.title.y=element_blank(),  #remove y axis labels,
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()  #remove y axis ticks
  ) +
  labs(title = "genes per cell")

percent_mito_per_cell <- ggplot(study_cell_stats, aes(x = `percent.mt`, y = sample_id, fill = study, group = sample_id)) +
  geom_density_ridges() +
  theme(
    axis.title.y=element_blank(),  #remove y axis labels,
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()  #remove y axis ticks
  ) +
  labs(title = "percent mito per cell") +
  xlim(0, 25)

library(patchwork)

study_stats_patchwork <- umis_per_cell + genes_per_cell + percent_mito_per_cell + plot_layout(guides = 'collect')

ggsave("results/study_stats.pdf", study_stats_patchwork)

# plot cell numbers ------------------------------

all_study_metadata <- read_csv("~/single_cell_projects/resources/external_rb_scrnaseq_proj/results/metadata_by_study.csv")

cells_per_sample <-
  all_study_metadata %>%
  group_by(study) %>%
  dplyr::count(sample_id) %>%
  identity()

write_csv(cells_per_sample, "results/cells_per_sample.csv")




