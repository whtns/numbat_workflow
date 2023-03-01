#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
library(cellpypes)
library(glue)
library(Seurat)
library(seuratTools)

celltype_markers <- read_csv("doc/celltype_marker_genes.csv")

# functions ------------------------------

plot_phylo_w_celltypes <- function(nb, myseu, myannot, mytitle, ...) {
  # browser()
  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    identity()

  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")

  mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = "GT_opt",
    ...
  ) +
    labs(title = mytitle)
}

safe_plot_phylo <- safely(plot_phylo_w_celltypes, otherwise = NA_real_)

compare_cell_vectors <- function(seu, cells.1 = NULL, cells.2 = NULL){
  mytibble <- list("group_1" = cells.1,
                   "group_2" = cells.2) %>%
    tibble::enframe("group", "cell") %>%
    tidyr::unnest(cell) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, mytibble)

  markers <- Seurat::FindMarkers(seu, ident.1 = "group_1", ident.2 = "group_2", group.by = "group")

  return(markers)

}

plot_subset_numbat <- function(nb, myannot, retained_cells, ...){
  nb2 <- nb$clone()

  nb2$joint_post <- nb2$joint_post[nb2$joint_post$cell %in% retained_cells,]

  nb2$plot_phylo_heatmap(
    clone_bar = FALSE,
    p_min = 0.2,
    annot = myannot,
    sort_by = "GT_opt",
    show_phylo = FALSE,
    ...
  )
}

# run_gsea(bg.genes = bg.genes, stats = stats.1,
#          category = 'C5', subcategory = 'BP',
#          out.dir = outdir, plot.title = 'GO',
#          file.prefix = file.prefix.1, n = 30)

suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(dplyr))

plot_go <-  function(gsea.res, n = 20, ...){
  # browser()
  gsea.res.order = order(gsea.res $padj, decreasing = FALSE)

  gsea.res$neg_log_padj = -log(gsea.res$padj)

  suppressPackageStartupMessages(library(ggplot2))

  # Format GOs
  # plot.table = head(gsea.res[gsea.res.order,], n = n)

  plot.table <-
    gsea.res %>%
    group_by(source) %>%
    dplyr::slice_head(n = n) %>%
    group_by(pathway) %>%
    dplyr::mutate(mean_NES = mean(NES)) %>%
    dplyr::arrange(mean_NES) %>%
    dplyr::mutate(pathway = str_remove(pathway, "GOBP_")) %>%
    dplyr::mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
    dplyr::mutate(pathway = factor(pathway, levels = unique(.[["pathway"]]))) %>%
    identity()

  # plot.table$pathway = sub('GOBP_', '', plot.table$pathway)
  # plot.table$pathway = gsub('_', ' ', plot.table$pathway)

  p = ggplot(plot.table,
             aes(x = NES, y = pathway)) +
    geom_point(aes(size = neg_log_padj, ...)) +
    theme_bw(base_size = 8) +
    ylab(NULL) +
    # scale_colour_gradient2(low = 'red',
    #                        mid = 'lightgrey',
    #                        high = 'blue',
    #                        midpoint = 0.05,
    #                        limits = c(0,0.1),
    #                        oob = scales::squish) +
    NULL

  return(p)
}

run_gsea <- function(bg.genes, stats, category = "C5", subcategory = "BP", n = 30){

  stats <- stats[c("symbol", "avg_log2FC")] %>%
    tibble::deframe()

  # Fetch geneset
  geneSets = msigdbr::msigdbr(species = 'Homo sapiens', category = category, subcategory = subcategory)
  geneSets = geneSets[geneSets$human_gene_symbol %in% bg.genes,]
  m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)

  # Run GSEA
  gsea = fgsea(pathways = m_list, stats = stats, minSize = 10, eps = 0.0)

  return(gsea)
}

process_diffex <- function(df, leukocyte_genes = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DPB1", "LST1", "AIF1", "PSMB8", "HLA-DQB1")){

  mycc_genes <-
    Seurat::cc.genes %>%
    tibble::enframe("phase", "symbol") %>%
    unnest(symbol) %>%
    identity()

  df %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(mycc_genes, by = "symbol") %>%
    dplyr::left_join(annotables::grch38[c("symbol", "chr", "description")], by = "symbol") %>%
    dplyr::filter(!symbol %in% leukocyte_genes) %>%
    dplyr::filter(!str_detect(chr, "CHR_"))
}

# preamble end ------------------------------

pdf("results/annotation_attempt.pdf")


# SRR13633759 normal control ------------------------------
study = "collin"
sample_id = "SRR13633759"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))

saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

# myannot = nb$clone_post[,c("cell", "GT_opt")]
#
# safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE)


# SRR13633760 normal control ------------------------------
study = "collin"
sample_id = "SRR13633760"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))

saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

# myannot = nb$clone_post[,c("cell", "GT_opt")]
#
# safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE)

# SRR13633761 FAILED ------------------------------
study = "collin"
sample_id = "SRR13633761"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE)


# SRR13633762 FAILED------------------------------
study = "collin"
sample_id = "SRR13633762"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE)

# # SRR13884244 FAILED------------------------------
# study = "wu"
# # sample_id = "SRR13884244"
# # seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
# #
# # # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# # #
# # # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
# #
# # nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
# #
# #
# #
# # myannot = nb$clone_post[,c("cell", "GT_opt")]
# #
# #
# #   safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#
# # SRR13884246 FAILED ------------------------------
# # too many SCNAs
#
# study = "wu"
# sample_id = "SRR13884246"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
#
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# myannot = nb$clone_post[,c("cell", "GT_opt")]
#
# safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#
# # SRR13884247 focal RB1-, focal MYCNA ------------------------------
# # focal RB deletion?
# # MYCNA?
#
# study = "wu"
# sample_id = "SRR13884247"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
#
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# myannot = nb$clone_post[,c("cell", "GT_opt")]
#
# safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#
# # SRR13884248 FAILED ------------------------------
# # focal deletion of 13q; no other SCNAs
#
# study = "wu"
# sample_id = "SRR13884248"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
#
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# myannot = nb$clone_post[,c("cell", "GT_opt")]
#
# safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)

# # SRR14800536 FAILED ------------------------------
# # no subclonal SCNAs
#
# sample_id = "SRR14800536"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
# #
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# # #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
#
#
# myannot = nb$clone_post[,c("cell", "GT_opt")]
#
#
#   safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#
# # SRR14800539 FAILED ------------------------------
# # no subclonal SCNAs
#
# sample_id = "SRR14800539"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
#
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# # #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]
#
#
#   safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#
#
# # SRR14800541 FAILED ------------------------------
# # too many spurious SCNAs
# # likely very few normal cells
#
# sample_id = "SRR14800541"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
#
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# # #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]
#
# safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#
# seu <- Seurat::AddMetaData(seu, tibble::deframe(nb$clone_post[,c("cell", "clone_opt")]), col.name = "clone_opt")
#
# DimPlot(seu, group.by = "type") +    labs(title = sample_id)
#
# seu <- seu[,!seu$type %in% c("Red Blood Cells")]
#
# DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)
#
# DimPlot(seu, group.by = "Phase") + labs(title = sample_id)
#
# # clone 1 vs. clone 2
# # negative means lower in clone 1; positive higher in clone 1
#
# clone_markers <- FindMarkers(seu, ident.1 = "2", group.by = "clone_opt") %>%
#   process_diffex()
#
#
# # SRR17960482 FAILED  ------------------------------
# # very weird looking
#
# sample_id = "SRR17960482"
# seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
#
# # nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# # #
# # saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))
#
# myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]
#
#
#   safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)
#

# 1q+ start ------------------------------

#* SRR13884243 wu 1q+ ------------------------------
# small quantity of 1qNULL cells

study = "wu"
sample_id = "SRR13884243"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884243_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

# seu <- seu[,!seu$type %in% c("Red Blood Cells")]

DimPlot(seu, group.by = "clone_opt")

DimPlot(seu, group.by = "gene_snn_res.0.2")

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)


# clone 1 vs. clone 2
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR13884243_one_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2", "3"), group.by = "clone_opt") %>%
  process_diffex()

SRR13884243_one_q_markers_gsea <- run_gsea(rownames(seu), SRR13884243_one_q_markers)

SRR13884243_one_q_markers_go_plot <- plot_go(SRR13884243_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)


marker_plot <- plot_markers(seu, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

plot_markers(seu, metavar = "gene_snn_res.0.2", marker_method = "genesorteR", return_plotly = FALSE)


obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_gene
)

pype <- obj %>%
  rule("macrophage1", "HLA-C", ">", 1) %>%
  # rule("macrophage2", "HLA-B", ">", 1, parent =  "macrophage1") %>%
  # rule("macrophage3", "HLA-C", ">", 1, parent =  "macrophage1") %>%
  # rule("macrophage3", "HLA-DRA", ">", 1, parent =  "macrophage1") %>%
  rule("macrophage3", "HLA-DPB1", ">", 1, parent =  "macrophage1") %>%
  # rule("macrophage4", "LST1", ">", 1, parent =  "macrophage1") %>%
  # plot_last() %>%
  identity()

plot_classes(pype)+ggtitle("SRR14800543 annotated with cellpypes")

SRR13884243_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")



matching_genes <- intersect(SRR14800534_one_q_markers$symbol,
                            SRR13884243_one_q_markers$symbol)

gsea_1q <-
  list(
    "SRR13884242" = SRR13884242_one_q_markers_gsea,
    "SRR13884243" = SRR13884243_one_q_markers_gsea,
    "SRR14800534" = SRR14800534_one_q_markers_gsea,
    "SRR14800537" = SRR14800537_one_q_markers_gsea,
    "SRR17960481" = SRR17960481_one_q_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  identity()

#* SRR14800534 !!! 1q+------------------------------
# large fraction of cells lacking 1q+
# almost all of the top genes are g2m genes down in clone 1!!
# 1q+ correlated with G2M?
# is this the same sample as SRR13884243?

study = "yang"
sample_id = "SRR14800534"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.95)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  # dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR14800534_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

# clone 1 vs. clone 2
# negative means lower in clone 1; positive higher in clone 1

# clone 1 vs. all
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR14800534_one_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = "2", group.by = "clone_opt") %>%
  process_diffex()

SRR14800534_one_q_markers_gsea <- run_gsea(rownames(seu), SRR14800534_one_q_markers)

SRR14800534_one_q_markers_go_plot <- plot_go(SRR14800534_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR14800534_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* SRR17960481 !!! 1q+, 6p+ ------------------------------
# notable fraction of 1q NULL cell
# PLK1 negative in clone 2 lacking 1q+

study = "field"
sample_id = "SRR17960481"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR17960481_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1", "2") ~ "scna_null",
    clone_opt %in% c("3") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

SRR17960481_6p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)

DimPlot(seu, group.by = "Phase") + labs(title = sample_id)


# clone 1,2 vs. 3
# clones 1,2 lack 1q+
# negative means higher in clone 3 (1q+); positive higher in clone 2 (1qNULL)
SRR17960481_one_q_markers <- FindMarkers(seu, ident.1 = c("2"), ident.2 = "3", group.by = "clone_opt") %>%
  process_diffex()

SRR17960481_one_q_markers_gsea <- run_gsea(rownames(seu), SRR17960481_one_q_markers)

SRR17960481_one_q_markers_go_plot <- plot_go(SRR17960481_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)


# clone 1 vs. 2,3
# clone 1 lacks 6p+
# negative means lower in (6pNULL); positive lower in (6p+)
SRR17960481_six_p_markers <- FindMarkers(seu, ident.1 = c("1"), ident.2 = c("2", "3"), group.by = "clone_opt") %>%
  process_diffex()

SRR17960481_six_p_markers_gsea <- run_gsea(rownames(seu), SRR17960481_six_p_markers)

SRR17960481_six_p_markers_go_plot <- plot_go(SRR17960481_six_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR17960481_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")


gsea_6p <- list(
  "SRR17960481" = SRR17960481_six_p_markers_gsea,
  "SRR17960484" = SRR17960484_six_p_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(pathway) %>%
  identity()

#* SRR14800537 !!! 1q+, 16q-------------------------------
# large fraction of cells missing 1q+/16q-
# few cell cycle related genes in diffex of clone 2 (1qNULL/16qNULL)
# lots of RB marker genes

study = "yang"
sample_id = "SRR14800537"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]

safe_plot_phylo(nb,
                seu,
                myannot,
                sample_id,
                clone_bar = FALSE,
                p_min = 0.7)

# 1q+ diffex subset nb phylo plot

clone_3_1q_null_cells <-
  nb$clone_post %>%
  dplyr::full_join(nb$joint_post, by = c("cell")) %>%
  dplyr::filter(clone_opt == "3") %>%
  dplyr::filter(CHROM == "1") %>%
  dplyr::filter(cnv_state == "loh", p_cnv.y < 0.1) %>%
  dplyr::pull(cell) %>%
  identity()

clone_4_cells <-
  nb$clone_post %>%
  dplyr::filter(clone_opt == "4") %>%
  dplyr::pull(cell) %>%
  identity()

clone_3_1q_clone_4_phylo_plot <- plot_subset_numbat(nb, myannot, c(clone_3_1q_null_cells, clone_4_cells))

# 16q- diffex subset nb phylo plot

clone_2_16q_null_cells <-
  nb$clone_post %>%
  dplyr::full_join(nb$joint_post, by = c("cell")) %>%
  dplyr::filter(clone_opt == "2") %>%
  dplyr::filter(CHROM == "16") %>%
  dplyr::filter(cnv_state == "del", p_cnv.y < 0.1) %>%
  dplyr::pull(cell) %>%
  identity()

clone_4_cells <-
  nb$clone_post %>%
  dplyr::filter(clone_opt == "4") %>%
  dplyr::pull(cell) %>%
  identity()

plot_subset_numbat(nb, myannot, c(clone_2_16q_null_cells, clone_4_cells))

SRR14800537_16q_diffex <-
  compare_cell_vectors(seu,
                       cells.1 = clone_2_16q_null_cells,
                       cells.2 = clone_4_cells)


SRR14800537_1q_diffex <-
  compare_cell_vectors(seu,
                       cells.1 = clone_3_1q_null_cells,
                       cells.2 = clone_4_cells)

oneq_null_v_sixteenq_null_diffex <-
  compare_cell_vectors(seu,
                       cells.1 = clone_2_16q_null_cells,
                       cells.2 = clone_3_1q_null_cells
  )


nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  # dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR14800537_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("3") ~ "scna_null",
    clone_opt %in% c("1", "2", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

SRR14800537_16q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("2") ~ "scna_null",
    clone_opt %in% c("1", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

seu <- seu[,!seu$type %in% c("Red Blood Cells")]

DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)

DimPlot(seu, group.by = "Phase") + labs(title = sample_id)

# clone 3 vs. clone 4
# negative means lower in clone 3; positive higher in clone 3 (1qNULL)
# clone 3 defferentiated relative to clone 4; lower expresison of cone pr genes
SRR14800537_one_q_markers <- FindMarkers(seu, ident.1 = "3", ident.2 = "4", group.by = "clone_opt") %>%
  process_diffex()

SRR14800537_one_q_markers_gsea <- run_gsea(rownames(seu), SRR14800537_one_q_markers)

SRR14800537_one_q_markers_go_plot <- plot_go(SRR14800537_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

# clone 2 vs. clone 4
# negative means lower in clone 2; positive higher in clone 2 (16qNULL)
# clone 2 has increased expression of TFF1
SRR14800537_sixteen_q_markers <- FindMarkers(seu, ident.1 = "2", ident.2 = "4", group.by = "clone_opt") %>%
  process_diffex()

SRR14800537_sixteen_q_markers_gsea <- run_gsea(rownames(seu), SRR14800537_sixteen_q_markers)

SRR14800537_sixteen_q_markers_go_plot <- plot_go(SRR14800537_sixteen_q_markers_gsea,
                                                 n = 20) +
  labs(title = sample_id)

SRR14800537_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* SRR13884242 wu 1q+, 16q-------------------------------
study = "wu"
sample_id = "SRR13884242"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt", "clone_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.99)

clone_1_1q_null_cells <-
  nb$clone_post %>%
  dplyr::full_join(nb$joint_post, by = c("cell")) %>%
  dplyr::filter(clone_opt == "1") %>%
  dplyr::filter(CHROM == "1") %>%
  dplyr::filter(cnv_state == "amp", p_cnv.y < 0.1) %>%
  dplyr::pull(cell) %>%
  identity()

clone_234_cells <-
  nb$clone_post %>%
  dplyr::filter(clone_opt %in% c("2", "3", "4")) %>%
  dplyr::pull(cell) %>%
  identity()

clone_1_1q_clone_234_phylo_plot <- plot_subset_numbat(nb, myannot, c(clone_1_1q_null_cells, clone_234_cells)) +
  labs(title = sample_id)

down_pdf(clone_1_1q_clone_234_phylo_plot, 8, 8)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884242_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

# clone 1 vs. all
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR13884242_one_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2", "3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR13884242_one_q_markers_gsea <- run_gsea(rownames(seu), SRR13884242_one_q_markers)

SRR13884242_one_q_markers_go_plot <- plot_go(SRR13884242_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

# clone 1 vs. all
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR13884242_sixteen_q_markers <- SRR13884242_one_q_markers

SRR13884242_sixteen_q_markers_gsea <- run_gsea(rownames(seu), SRR13884242_sixteen_q_markers)

SRR13884242_sixteen_q_markers_go_plot <- plot_go(SRR13884242_sixteen_q_markers_gsea,
                                                 n = 20) +
  labs(title = sample_id)

SRR13884242_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* SRR17960484 !!! 6p+, 1q+ ------------------------------
# very few 1qNULL clone 1 diffex 1 v. 2
# possibly some 6pNULL clone 1/2 diffex 1/2 v. 3/4

sample_id = "SRR17960484"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR17960484_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

SRR17960484_6p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("2") ~ "scna_null",
    clone_opt %in% c("1", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

seu <- seu[,!seu$type %in% c("Red Blood Cells")]

DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)

DimPlot(seu, group.by = "Phase") + labs(title = sample_id)

# clone 1 vs. clone 2
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR17960484_one_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2", "3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR17960484_one_q_markers_gsea <- run_gsea(rownames(seu), SRR17960484_one_q_markers)

SRR17960484_one_q_markers_go_plot <- plot_go(SRR17960484_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

# clone 2 vs. clone 3/4
# clone 2 lacks 16q-
# negative means lower in clone 1/2; positive higher in clone 3/4
SRR17960484_six_p_markers <- FindMarkers(seu, ident.1 = c("2"), ident.2 = c("3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR17960484_six_p_markers_gsea <- run_gsea(rownames(seu), SRR17960484_six_p_markers)

SRR17960484_six_p_markers_go_plot <- plot_go(SRR17960484_six_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR17960484_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* cohort ------------------------------

one_q_gseas <- list(SRR13884242 = SRR13884242_one_q_markers_gsea,
                    SRR13884243 = SRR13884243_one_q_markers_gsea,
                    SRR14800534 = SRR14800534_one_q_markers_gsea,
                    SRR14800537 = SRR14800537_one_q_markers_gsea,
                    SRR17960481 = SRR17960481_one_q_markers_gsea,
                    SRR17960484 = SRR17960484_one_q_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  identity()

one_q_go_plot <- plot_go(one_q_gseas,
                         n = 20, color = source, shape = source) +
  labs(title = "1q+ enrichment") +
  scale_y_discrete(     labels=function(x)str_wrap(x, width = 60))


# 1q+ end ------------------------------

# 2p+ start ------------------------------

# SRR13884240 wu 2p+------------------------------
study = "wu"
sample_id = "SRR13884240"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884240_2p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

# clone 1 vs. clone 2,3,4
# negative means lower in clone 1; positive higher in clone 2 (16qNULL)
SRR13884240_two_p_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2", "3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR13884240_two_p_markers_gsea <- run_gsea(rownames(seu), SRR13884240_two_p_markers)

SRR13884240_two_p_markers_go_plot <- plot_go(SRR13884240_two_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR13884240_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

# SRR13884241 wu 2p+ ------------------------------
study = "wu"
sample_id = "SRR13884241"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884241_2p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

# clone 1 vs. clone 2
# negative means lower in clone 1; positive higher in clone 2 (16qNULL)
SRR13884241_two_p_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2", "3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR13884241_two_p_markers_gsea <- run_gsea(rownames(seu), SRR13884241_two_p_markers)

SRR13884241_two_p_markers_go_plot <- plot_go(SRR13884241_two_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR13884241_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")


gsea_2p <-
  list(
    "SRR13884240" = SRR13884240_two_p_markers_gsea,
    "SRR13884241" = SRR13884241_two_p_markers_gsea,
    "SRR13884249" = SRR13884249_two_p_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(pathway) %>%
  identity()

# SRR13884249 wu 2p+ ------------------------------
# large fraction of cells lacking 2p+

study = "wu"
sample_id = "SRR13884249"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.5)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  # dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884249_2p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

# clone 1 vs. clone 2
# negative means lower in clone 1; positive higher in clone 2 (16qNULL)
SRR13884249_two_p_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2"), group.by = "clone_opt") %>%
  process_diffex()

SRR13884249_two_p_markers_gsea <- run_gsea(rownames(seu), SRR13884249_two_p_markers)

SRR13884249_two_p_markers_go_plot <- plot_go(SRR13884249_two_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR13884249_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* cohort ------------------------------

two_p_gseas <- list(
  SRR13884240 = SRR13884240_two_p_markers_gsea,
  SRR13884241 = SRR13884241_two_p_markers_gsea,
  SRR13884249 = SRR13884249_two_p_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  identity()

two_p_go_plot <- plot_go(two_p_gseas,
                         n = 20, color = source, shape = source) +
  labs(title = "2p+ enrichment") +
  scale_y_discrete(     labels=function(x)str_wrap(x, width = 60))

# 2p+ end ------------------------------

# 6p+ start ------------------------------

# SRR17960481 !!! 1q+, 6p+ ------------------------------
# notable fraction of 1q NULL cell
# PLK1 negative in clone 2 lacking 1q+

study = "field"
sample_id = "SRR17960481"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR17960481_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1", "2") ~ "scna_null",
    clone_opt %in% c("3") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

SRR17960481_6p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)

DimPlot(seu, group.by = "Phase") + labs(title = sample_id)


# clone 1 vs. 2,3
# clone 1 lacks 6p+
# negative means lower in (6pNULL); positive lower in (6p+)
SRR17960481_six_p_markers <- FindMarkers(seu, ident.1 = c("1"), ident.2 = c("2", "3"), group.by = "clone_opt") %>%
  process_diffex()

SRR17960481_six_p_markers_gsea <- run_gsea(rownames(seu), SRR17960481_six_p_markers)

SRR17960481_six_p_markers_go_plot <- plot_go(SRR17960481_six_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR17960481_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")


gsea_6p <- list(
  "SRR17960481" = SRR17960481_six_p_markers_gsea,
  "SRR17960484" = SRR17960484_six_p_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(pathway) %>%
  identity()

# SRR17960484 !!! 6p+, 1q+ ------------------------------
# very few 1qNULL clone 1 diffex 1 v. 2
# possibly some 6pNULL clone 1/2 diffex 1/2 v. 3/4

sample_id = "SRR17960484"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR17960484_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

SRR17960484_6p_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("2") ~ "scna_null",
    clone_opt %in% c("1", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

seu <- seu[,!seu$type %in% c("Red Blood Cells")]

DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)

DimPlot(seu, group.by = "Phase") + labs(title = sample_id)


# clone 2 vs. clone 3/4
# clone 2 lacks 16q-
# negative means lower in clone 1/2; positive higher in clone 3/4
SRR17960484_six_p_markers <- FindMarkers(seu, ident.1 = c("2"), ident.2 = c("3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR17960484_six_p_markers_gsea <- run_gsea(rownames(seu), SRR17960484_six_p_markers)

SRR17960484_six_p_markers_go_plot <- plot_go(SRR17960484_six_p_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

SRR17960484_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* cohort ------------------------------

six_p_gseas <- list(
  SRR17960481 = SRR17960481_six_p_markers_gsea,
  SRR17960484 = SRR17960484_six_p_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  identity()

six_p_go_plot <- plot_go(six_p_gseas,
                         n = 20, color = source, shape = source) +
  labs(title = "6p+ enrichment")

# 16q- start ------------------------------

# SRR14800535 !!! 16q-------------------------------
# very small fraction of cells with 16q-
# PLK1 diff expressed between clones; lower in clone 1

study = "yang"
sample_id = "SRR14800535"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  # dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR14800535_16q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

seu <- seu[,!seu$type %in% c("Red Blood Cells")]

DimPlot(seu, group.by = "clone_opt")

DimPlot(seu, group.by = "gene_snn_res.0.2")

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)

# clone 1 vs. clone 2
# clone 1 lacks 16q-
# negative means lower in clone 1; positive higher in clone 1
SRR14800535_sixteen_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2"), group.by = "clone_opt") %>%
  process_diffex()

SRR14800535_sixteen_q_markers_gsea <- run_gsea(rownames(seu), SRR14800535_sixteen_q_markers)

SRR14800535_sixteen_q_markers_go_plot <- plot_go(SRR14800535_sixteen_q_markers_gsea,
                                                 n = 20) +
  labs(title = sample_id)

SRR14800535_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")



# SRR14800537 !!! 1q+, 16q-------------------------------
# large fraction of cells missing 1q+/16q-
# few cell cycle related genes in diffex of clone 2 (1qNULL/16qNULL)
# lots of RB marker genes

study = "yang"
sample_id = "SRR14800537"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
# #
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "clone_opt", "GT_opt")]

safe_plot_phylo(nb,
                seu,
                myannot,
                sample_id,
                clone_bar = FALSE,
                p_min = 0.7)

# 1q+ diffex subset nb phylo plot

clone_3_1q_null_cells <-
  nb$clone_post %>%
  dplyr::full_join(nb$joint_post, by = c("cell")) %>%
  dplyr::filter(clone_opt == "3") %>%
  dplyr::filter(CHROM == "1") %>%
  dplyr::filter(cnv_state == "loh", p_cnv.y < 0.1) %>%
  dplyr::pull(cell) %>%
  identity()

clone_4_cells <-
  nb$clone_post %>%
  dplyr::filter(clone_opt == "4") %>%
  dplyr::pull(cell) %>%
  identity()

clone_3_1q_clone_4_phylo_plot <- plot_subset_numbat(nb, myannot, c(clone_3_1q_null_cells, clone_4_cells))

# 16q- diffex subset nb phylo plot

clone_2_16q_null_cells <-
  nb$clone_post %>%
  dplyr::full_join(nb$joint_post, by = c("cell")) %>%
  dplyr::filter(clone_opt == "2") %>%
  dplyr::filter(CHROM == "16") %>%
  dplyr::filter(cnv_state == "del", p_cnv.y < 0.1) %>%
  dplyr::pull(cell) %>%
  identity()

clone_4_cells <-
  nb$clone_post %>%
  dplyr::filter(clone_opt == "4") %>%
  dplyr::pull(cell) %>%
  identity()

plot_subset_numbat(nb, myannot, c(clone_2_16q_null_cells, clone_4_cells))

SRR14800537_16q_diffex <-
  compare_cell_vectors(seu,
                       cells.1 = clone_2_16q_null_cells,
                       cells.2 = clone_4_cells)


SRR14800537_1q_diffex <-
  compare_cell_vectors(seu,
                       cells.1 = clone_3_1q_null_cells,
                       cells.2 = clone_4_cells)

oneq_null_v_sixteenq_null_diffex <-
  compare_cell_vectors(seu,
                       cells.1 = clone_2_16q_null_cells,
                       cells.2 = clone_3_1q_null_cells
  )


nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  # dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR14800537_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("3") ~ "scna_null",
    clone_opt %in% c("1", "2", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

SRR14800537_16q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("2") ~ "scna_null",
    clone_opt %in% c("1", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

seu <- seu[,!seu$type %in% c("Red Blood Cells")]

DimPlot(seu, group.by = "clone_opt") + labs(title = sample_id)

DimPlot(seu, group.by = "Phase") + labs(title = sample_id)

# clone 3 vs. clone 4
# negative means lower in clone 3; positive higher in clone 3 (1qNULL)
# clone 3 defferentiated relative to clone 4; lower expresison of cone pr genes
SRR14800537_one_q_markers <- FindMarkers(seu, ident.1 = "3", ident.2 = "4", group.by = "clone_opt") %>%
  process_diffex()

SRR14800537_one_q_markers_gsea <- run_gsea(rownames(seu), SRR14800537_one_q_markers)

SRR14800537_one_q_markers_go_plot <- plot_go(SRR14800537_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

# clone 2 vs. clone 4
# negative means lower in clone 2; positive higher in clone 2 (16qNULL)
# clone 2 has increased expression of TFF1
SRR14800537_sixteen_q_markers <- FindMarkers(seu, ident.1 = "2", ident.2 = "4", group.by = "clone_opt") %>%
  process_diffex()

SRR14800537_sixteen_q_markers_gsea <- run_gsea(rownames(seu), SRR14800537_sixteen_q_markers)

SRR14800537_sixteen_q_markers_go_plot <- plot_go(SRR14800537_sixteen_q_markers_gsea,
                                                 n = 20) +
  labs(title = sample_id)

SRR14800537_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

# SRR13884242 wu 1q+, 16q-------------------------------
study = "wu"
sample_id = "SRR13884242"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

myannot = nb$clone_post[,c("cell", "GT_opt", "clone_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.99)

clone_1_1q_null_cells <-
  nb$clone_post %>%
  dplyr::full_join(nb$joint_post, by = c("cell")) %>%
  dplyr::filter(clone_opt == "1") %>%
  dplyr::filter(CHROM == "1") %>%
  dplyr::filter(cnv_state == "amp", p_cnv.y < 0.1) %>%
  dplyr::pull(cell) %>%
  identity()

clone_234_cells <-
  nb$clone_post %>%
  dplyr::filter(clone_opt %in% c("2", "3", "4")) %>%
  dplyr::pull(cell) %>%
  identity()

clone_1_1q_clone_234_phylo_plot <- plot_subset_numbat(nb, myannot, c(clone_1_1q_null_cells, clone_234_cells)) +
  labs(title = sample_id)

down_pdf(clone_1_1q_clone_234_phylo_plot, 8, 8)

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884242_1q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2", "3", "4") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

# clone 1 vs. all
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR13884242_one_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2", "3", "4"), group.by = "clone_opt") %>%
  process_diffex()

SRR13884242_one_q_markers_gsea <- run_gsea(rownames(seu), SRR13884242_one_q_markers)

SRR13884242_one_q_markers_go_plot <- plot_go(SRR13884242_one_q_markers_gsea,
                                             n = 20) +
  labs(title = sample_id)

# clone 1 vs. all
# clone 1 lacks 1q+
# negative means lower in clone 1; positive higher in clone 1
SRR13884242_sixteen_q_markers <- SRR13884242_one_q_markers

SRR13884242_sixteen_q_markers_gsea <- run_gsea(rownames(seu), SRR13884242_sixteen_q_markers)

SRR13884242_sixteen_q_markers_go_plot <- plot_go(SRR13884242_sixteen_q_markers_gsea,
                                                 n = 20) +
  labs(title = sample_id)

SRR13884242_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

# SRR13884245 wu 16q- ------------------------------
# a very small number of cells with 16q-

study = "wu"
sample_id = "SRR13884245"
seu <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))

# nb = Numbat$new(out_dir = glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}"))
#
# saveRDS(nb, glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- nb$clone_post[,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

seu <- Seurat::AddMetaData(seu, nb_meta)

SRR13884245_16q_clone_by_cluster <- FetchData(seu, c("clone_opt", "gene_snn_res.0.2")) %>%
  dplyr::filter(!clone_opt == "") %>%
  dplyr::mutate(diffex_group = dplyr::case_when(
    clone_opt %in% c("1") ~ "scna_null",
    clone_opt %in% c("2") ~ "scna"
  )) %>%
  janitor::tabyl(`gene_snn_res.0.2`, diffex_group) %>%
  identity()

GT_by_cluster_plot <- ggplot(seu@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
  geom_bar(position="fill") +
  coord_flip()

cluster_plot <- DimPlot(seu, group.by = "gene_snn_res.0.2")
GT_plot <- DimPlot(seu, group.by = "GT_opt")

GT_by_cluster_plot / (cluster_plot + GT_plot)

DimPlot(seu, group.by = "type") +    labs(title = sample_id)

seu <- seu[,!seu$type %in% c("Microglia")]

DimPlot(seu, group.by = "clone_opt")

DimPlot(seu, group.by = "gene_snn_res.0.2")

myannot = nb$clone_post[,c("cell", "GT_opt")]

safe_plot_phylo(nb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.2)

# clone 1 vs. clone 2
# clone 1 lacks 16q-
# negative means lower in clone 1; positive higher in clone 1
SRR13884245_sixteen_q_markers <- FindMarkers(seu, ident.1 = "1", ident.2 = c("2"), group.by = "clone_opt") %>%
  process_diffex()


SRR13884245_sixteen_q_markers_gsea <- run_gsea(rownames(seu), SRR13884245_sixteen_q_markers)

SRR13884245_sixteen_q_markers_go_plot <- plot_go(SRR13884245_sixteen_q_markers_gsea,
                                                 n = 20) +
  labs(title = sample_id)

SRR13884245_dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2", split.by = "clone_opt")

#* cohort ------------------------------

sixteen_q_gseas <- list(
  "SRR14800535" = SRR14800535_sixteen_q_markers_gsea,
  "SRR14800537" = SRR14800537_sixteen_q_markers_gsea,
  "SRR13884242" = SRR13884242_sixteen_q_markers_gsea,
  "SRR13884245" = SRR13884245_sixteen_q_markers_gsea) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(pathway) %>%
  identity()
#
# matching_genes <- intersect(SRR14800535_sixteen_q_markers$symbol,
#                             SRR13884245_sixteen_q_markers$symbol)
#

sixteen_q_go_plot <- plot_go(sixteen_q_gseas,
                         n = 10, color = source, shape = source) +
  labs(title = "16q- enrichment")


# 16q- end ------------------------------

dev.off()

gseas <-
  list(
    one_q = one_q_gseas,
    two_p = two_p_gseas,
    six_p = six_p_gseas,
    sixteen_q = sixteen_q_gseas
    ) %>%
  dplyr::bind_rows(.id = "region") %>%
  group_by(pathway) %>%
  dplyr::mutate(num_regions = n_distinct(region)) %>%
  dplyr::arrange(desc(num_regions)) %>%
  identity()

all_go_plot <-
  plot_go(gseas, color = source)

positive_gseas <-
  gseas %>%
  dplyr::filter(NES > 0)

positive_go_plot <-
  plot_go(positive_gseas, color = source) +
  scale_y_discrete(     labels=function(x)str_wrap(x, width = 60))

negative_gseas <-
  gseas %>%
  dplyr::filter(NES < 0)

negative_go_plot <-
  plot_go(negative_gseas, color = source) +
  scale_y_discrete(     labels=function(x)str_wrap(x, width = 60))

# clone by clusters ------------------------------

clone_by_clusters <- ls(pattern = "*_clone_by_cluster") %>%
  set_names(.) %>%
  map(get) %>%
  dplyr::bind_rows(.id = "source") %>%
  # dplyr::filter(padj < 0.05) %>%
  tidyr::pivot_longer(-c("source", "gene_snn_res.0.2"), names_to = "scna", values_to = "num_cells") %>%
  dplyr::mutate(cluster = `gene_snn_res.0.2`) %>%
  dplyr::mutate(region = str_extract(source, "[0-9]*[p,q]")) %>%
  identity()

cluster_plot <- ggplot(clone_by_clusters, aes(cluster, num_cells, fill = scna)) +
  geom_bar(position = "fill", stat = "identity") +
  # scale_y_log10() +
  facet_wrap(~source) +
  NULL

# dimplots ------------------------------
dimplots <- ls(pattern = "*_dimplot") %>%
  set_names(.) %>%
  map(get) %>%
  identity()

# gseas ------------------------------
gsea_tables <- ls(pattern = "*_markers_gsea") %>%
  set_names(.) %>%
  map(get) %>%
  dplyr::bind_rows(.id = "source") %>%
  dplyr::filter(padj < 0.05) %>%
  identity()


# go plots ------------------------------

pdf("results/go_plots.pdf")
go_plots <- ls(pattern = "*_go_plot") %>%
  map(get)
go_plots
dev.off()



