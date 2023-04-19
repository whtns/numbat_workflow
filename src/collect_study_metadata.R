#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)
library(glue)
library(ggridges)


retrieve_cell_stats <- function(seu_path){

  seu <- readRDS(seu_path)

  stats = seu@meta.data[c("nCount_gene", "nFeature_gene", "percent.mt", "gene_snn_res.0.2")] %>%
    tibble::rownames_to_column("cell")

  return(stats)
}

seus <-
  dir_ls("output/seurat/", regexp = ".*[0-9]+_seu.rds") %>%
  set_names(str_extract(., "SRR[0-9]*"))


# collin ------------------------------

collin_cell_stats <- seus[c("SRR13633759", "SRR13633760", "SRR13633761", "SRR13633762")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# field ------------------------------

field_cell_stats <- seus[c("SRR17960480", "SRR17960481", "SRR17960482", "SRR17960483", "SRR17960484")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# wu ------------------------------

wu_cell_stats <- seus[c("SRR13884240", "SRR13884241", "SRR13884242", "SRR13884243", "SRR13884244", "SRR13884245", "SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# yang ------------------------------

yang_cell_stats <- seus[c("SRR14800534", "SRR14800535", "SRR14800536", "SRR14800537", "SRR14800538", "SRR14800539", "SRR14800540", "SRR14800541", "SRR14800542", "SRR14800543")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# combined ------------------------------

study_cell_stats <- dplyr::bind_rows(list("collin" = collin_cell_stats, "field" = field_cell_stats, "wu" = wu_cell_stats, "yang" = yang_cell_stats), .id = "study")

umis_per_cell <- ggplot(study_cell_stats, aes(x = nCount_gene, y = sample_id, fill = study, group = sample_id)) +
  geom_density_ridges() +
  scale_x_log10() +
  scale_y_discrete(limits=rev) +
  theme(
    axis.title.y=element_blank(),  #remove y axis labels,
    # axis.text.y=element_blank(),  #remove y axis labels
    # axis.ticks.y=element_blank()  #remove y axis ticks
  ) +
  labs(title = "UMIs per cell")

genes_per_cell <- ggplot(study_cell_stats, aes(x = nFeature_gene, y = sample_id, fill = study, group = sample_id)) +
  geom_density_ridges() +
  scale_x_log10() +
  scale_y_discrete(limits=rev) +
  theme(
    axis.title.y=element_blank(),  #remove y axis labels,
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()  #remove y axis ticks
  ) +
  labs(title = "genes per cell")

percent_mito_per_cell <- ggplot(study_cell_stats, aes(x = `percent.mt`, y = sample_id, fill = study, group = sample_id)) +
  geom_density_ridges() +
  scale_y_discrete(limits=rev) +
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




