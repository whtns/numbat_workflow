#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
library(cellpypes)

pdf("results/cell_pype_annotations.pdf")

# SRR13633759 ------------------------------
seu <- normal_seus$SRR13633759

obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_RNA
)

pype <- obj %>%
  rule("connective",           "COL1A1",   ">", 7)                    %>%
  rule("lymphocyte",           "HBB",   ">", 100)                    %>%
  rule("cornea",           "KRT12",   ">", 10)                    %>%
  rule("macrophage",           "RNASE1",   ">", 10)                    %>%
  rule("conjunctival epithelial",           "S100A9",   ">", 10)                    %>%
  rule("photoreceptor",           "AIPL1",   ">", 10)                    %>%
  rule("ECM",           "ANGPTL7",   ">", 10)                    %>%
  rule("plasma",           "HLA-B",   ">", 10)                    %>%
  rule("ganglion", "CRABP1", ">", 4) %>%
  # rule("Memory CD4+",  "S100A4", ">", 13,  parent="CD4+ T") %>%
  identity()

classes_plot <- plot_classes(pype)+ggtitle("SRR13633759 annotated with cellpypes")

celltypes <- classify(pype)

seu <- AddMetaData(seu, celltypes, col.name = "CellType_predict")

# SRR14800540 ------------------------------
seu <- normal_seus$SRR14800540

obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_gene
)

pype <- obj %>%
  rule("cone",           "ARR3",   ">", 1)                    %>%
  rule("rod",           "GNAT1",   ">", 6)                    %>%
  rule("muscle",           "RGS5",   ">", 4)                    %>%
  rule("lymphocyte",           "C1QB",   ">", 4)                    %>%
  # rule("plasma",           "HLA-B",   ">", 20)                    %>%
  # rule("Memory CD4+",  "S100A4", ">", 13,  parent="CD4+ T") %>%
  identity()

plot_classes(pype)+ggtitle("SRR14800540 annotated with cellpypes")


# SRR14800541 ------------------------------
seu <- normal_seus$SRR14800541

obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_gene
)

pype <- obj %>%
  rule("cone",           "ARR3",   ">", 1)                    %>%
  rule("rod",           "GNAT1",   ">", 6)                    %>%
  rule("muscle",           "RGS5",   ">", 4)                    %>%
  rule("lymphocyte",           "C1QB",   ">", 4)                    %>%
  # rule("plasma",           "HLA-B",   ">", 20)                    %>%
  # rule("Memory CD4+",  "S100A4", ">", 13,  parent="CD4+ T") %>%
  identity()

plot_classes(pype)+ggtitle("SRR14800541 annotated with cellpypes")


# SRR14800542 ------------------------------
seu <- normal_seus$SRR14800542

obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_gene
)

pype <- obj %>%
  rule("cone",           "ARR3",   ">", 1)                    %>%
  rule("rod",           "GNAT1",   ">", 6)                    %>%
  rule("muscle",           "RGS5",   ">", 4)                    %>%
  rule("lymphocyte",           "C1QB",   ">", 4)                    %>%
  rule("t/nk cells",           "GZMA",   ">", 4)                    %>%
  rule("macrophage",           "TYROBP",   ">", 4)                    %>%
  # rule("plasma",           "HLA-B",   ">", 20)                    %>%
  # rule("Memory CD4+",  "S100A4", ">", 13,  parent="CD4+ T") %>%
  identity()

plot_classes(pype)+ggtitle("SRR14800542 annotated with cellpypes")


# SRR14800543 ------------------------------
seu <- normal_seus$SRR14800543

obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_gene
)

pype <- obj %>%
  rule("cone",           "ARR3",   ">", 1)                    %>%
  rule("rod",           "GNAT1",   ">", 6)                    %>%
  rule("muscle",           "RGS5",   ">", 4)                    %>%
  rule("lymphocyte",           "C1QB",   ">", 4)                    %>%
  rule("t/nk cells",           "GZMA",   ">", 4)                    %>%
  rule("macrophage",           "TYROBP",   ">", 4)                    %>%
  rule("retinoblastoma",           "MYCN",   ">", 3)                    %>%
  # rule("plasma",           "HLA-B",   ">", 20)                    %>%
  # rule("Memory CD4+",  "S100A4", ">", 13,  parent="CD4+ T") %>%
  identity()

plot_classes(pype)+ggtitle("SRR14800543 annotated with cellpypes")

dev.off()

# merged ------------------------------
plot_patchwork <- function(seu){

  marker_plot <- plot_markers(seu, metavar = "gene_snn_res.0.05", marker_method = "presto", return_plotly = FALSE, num_markers = 10)

  umap_plot <- DimPlot(seu, group.by = "gene_snn_res.0.05")

  mypatchwork <- umap_plot + marker_plot + plot_layout(widths = c(3,1))

  return(mypatchwork)
}


test0 <- plot_patchwork(seu)

seu <- merged_normal_seu0

obj <- list(
  raw      =SeuratObject::GetAssayData(seu, "counts"),
  neighbors=as(seu@graphs[["gene_nn"]], "dgCMatrix")>.1, # sometims "wknn"
  embed    =FetchData(seu, c("UMAP_1", "UMAP_2")),
  totalUMI = seu$nCount_gene
)

pype <- obj %>%
  rule("neuronal", "CCL5", ">", 2) %>%
  rule("cone", "ARR3", ">", 2) %>%
  rule("rpe", "DCT", ">", 2) %>%
  rule("bipolar cell", "GNG13", ">", 2) %>%
  rule("mesenchymal", "APOA1", ">", 2) %>%
  rule("muller glia", "PTGDS", ">", 2) %>%
  rule("rod", "GNAT1", ">", 6) %>%
  rule("retinoblastoma", "HES6", ">", 2) %>%
  rule("cornea", "KRT5", ">", 3) %>%
  rule("macrophage", "HLA-DRA", ">", 2) %>%
  rule("epithelial", "FN1", ">", 2) %>%
  # plot_last() %>%
  identity()

plot_classes(pype)+ggtitle("SRR14800543 annotated with cellpypes")


merged_normal_seu <- readRDS("output/seurat/merged_normal_seu.rds")


