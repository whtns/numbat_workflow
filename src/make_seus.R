#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(glue)
library(patchwork)

library(Seurat)
library(seuratTools)
library(infercnv)

normal_reference_mat <-
  readRDS("output/infercnv/reference_counts.rds")
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
  RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

make_seus_from_cellranger <- function(sample_path, normal_seu) {
  # browser()
  seu_path <-
    path("output/seurat", paste0(path_file(sample_path), "_seu.rds"))

  print(seu_path)

  mypath <- fs::path(sample_path, "outs/filtered_feature_bc_matrix")

  count_mat <- Seurat::Read10X(mypath)

  seu <- Seurat::CreateSeuratObject(count_mat, assay = "gene") %>%
    RenameCells(new.names = str_replace(colnames(.), "-", "."))

  seu <-
    seuratTools::clustering_workflow(seu, resolution = c(0.2, 0.4))

  saveRDS(seu, seu_path)
  print(glue::glue("saved {seu_path}"))

  return(seu)

  # seu_merged <- merge(seu, normal_seu)
  #
  # seu_merged <- infercnv::add_to_seurat(seu_merged, fs::path("output/infercnv", path_file(sample_path)))
  #
  # seu_w_cnv <- seu_merged[,!grepl("normal", colnames(seu_merged))]
  #
  # seu_w_cnv <- seuratTools::clustering_workflow(seu_w_cnv, resolution = c(0.2, 0.4))
  #
  # seu_cnv_path <- path("output/seurat", paste0(path_file(sample_path), "_cnv_seu.rds"))
  #
  # saveRDS(seu_w_cnv, seu_cnv_path)
  #
  # seu_cnv_path

}

study_dirs <-
  c("collin_et_al", "yang_et_al", "wu_et_al", "field_et_al")

seu_paths <-
  purrr::map(paste0("output/seurat/", study_dirs),
             ~ fs::dir_ls(.x, glob = "*rds")) %>%
  unlist()

normal_seu_paths <-
  c(
    "output/seurat/collin_et_al/SRR13633759.rds",
    "output/seurat/yang_et_al/SRR14800540_seu.rds",
    "output/seurat/yang_et_al/SRR14800541_seu.rds",
    "output/seurat/yang_et_al/SRR14800542_seu.rds",
    "output/seurat/yang_et_al/SRR14800543_seu.rds"
  ) %>%
  purrr::set_names(stringr::str_extract(path_file(.), "SRR[0-9]*"))

normal_seus <- purrr::map(normal_seu_paths, readRDS)

test0 <- purrr::map(normal_seus, RenameAssays, gene = "RNA")

test1 <- purrr::map(test0, clustering_workflow, resolution = 0.2, assay = "RNA")

tumor_seu_paths <- seu_paths[!seu_paths %in% normal_seu_paths]

plot_patchwork <- function(seu_path){

    sample_id <- stringr::str_extract(path_file(seu_path), "SRR[0-9]*")

    seu <- readRDS(seu_path)

    marker_plot <- plot_markers(seu, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

    umap_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2") + labs(title = sample_id)

    mypatchwork <- umap_plot + marker_plot + plot_layout(widths = c(3,1))

    ggsave(glue::glue("results/{sample_id}_patchwork.pdf"), mypatchwork, width = 10, height = 8)

    return(mypatchwork)
}

purrr::map(tumor_seu_paths, plot_patchwork)

purrr::map(normal_seu_paths[4], plot_patchwork)

marker_plots <- purrr::map(normal_seus, ~plot_markers(.x, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE))

umap_plots <- purrr::imap(seus, ~(DimPlot(.x, group.by = "gene_snn_res.0.2") + labs(title = .y)))

library(patchwork)

patchworks <- purrr::map2(umap_plots, marker_plots, ~.x + .y + plot_layout(widths = c(3,1)))

pdf("results/patchworks.pdf", width = 8)
patchworks
dev.off()


infercnv_obj <- readRDS("output/infercnv/SRR13884240/run.final.infercnv_obj")

seu_40 <- readRDS("output/seurat/SRR13884240_seu.rds") %>%
	RenameCells(new.names = str_replace(colnames(.), "-", "."))

normal_reference_mat <- readRDS("output/infercnv/reference_counts.rds")
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
	RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

seu0 <- merge(seu_40, normal_seu)

seu_cells <- rownames(seu0@meta.data)

infercnv_cells <- colnames(infercnv_obj@expr.data)


# debug(infercnv::add_to_seurat)
seu1 <- infercnv::add_to_seurat(seu0, "output/infercnv/SRR13884240/")

seu2 <- seu1[,!grepl("normal", colnames(seu1))]

seu2 <- seurat_preprocess(seu2) %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)

has_cnv_cols <- str_subset(colnames(seu2@meta.data), "has_cnv*")

cnv_plots <- purrr::map(has_cnv_cols, ~DimPlot(seu2, group.by = .x))

pdf("~/tmp/cnvplots.pdf")
cnv_plots
dev.off()
