#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
library(Seurat)
library(seuratTools)
library(patchwork)
library(msigdbr)
library(glue)

hallmark_gene_sets = msigdbr(species = "Homo sapiens", category = "H")

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


# plot variability at SCNA
plot_variability_at_SCNA <- function(phylo_plot_output, chrom = "1", p_min = 0.9){
	test0 <-
		phylo_plot_output$data %>%
		dplyr::mutate(seg = factor(seg, levels = str_sort(unique(seg), numeric = TRUE))) %>%
		# dplyr::filter(CHROM == chrom) %>%
		identity()

	p_cnv_plot <- ggplot(test0, aes(x = cell_index, y = p_cnv, color = cnv_state, alpha = p_cnv)) +
		geom_point(size = 0.1) +
	  geom_hline(aes(yintercept = p_min)) +
	  scale_x_reverse() +
		facet_wrap(~seg)

	return(p_cnv_plot)
}

read_regress_save <- function(seu_path){
	seu <- readRDS(seu_path)

	seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seu))

	# cell cycle effects strongly mitigated in PCA
	seu <- RunPCA(seu, features = VariableFeatures(seu), nfeatures.print = 10)

	seu <- FindNeighbors(seu, dims = 1:10)
	seu <- FindClusters(seu, resolution = 0.15)

	seu <- RunUMAP(seu, dims = 1:10, min.dist = 0.01)

	saveRDS(seu, seu_path)
}

read_unregress_cc_save <- function(seu_path){
	seu <- readRDS(seu_path)

	seu <- seuratTools::clustering_workflow(seu, resolution = c(0.2, 0.4))

	# cell cycle effects strongly mitigated in PCA
	# seu <- seuratTools::seurat_reduce_dimensions(seu)

	# seu <- RunPCA(seu, features = VariableFeatures(seu), nfeatures.print = 10)
	#
	# seu <- RunUMAP(seu, dims = 1:30)

	# saveRDS(seu, seu_path)
}

annotate_seu_with_rb_subtype_gene_expression <- function(seu){
	subtype_genes <- list(
		gp1 = c("EGF", "TPBG", "GUCA1C", "GUCA1B", "GUCA1A", "GNAT2", "GNGT2", "ARR3", "PDE6C", "PDE6H", "OPN1SW"),
		gp2 = c("TFF1", "CD24", "EBF3", "GAP43", "STMN2", "POU4F2", "SOX11", "EBF1", "DCX", "ROBO1", "PCDHB10")
	)

	seu <- AddModuleScore(
		object = seu,
		features = subtype_genes,
		ctrl = 5,
		name = 'exprs_gp'
	)

	subtype_hallmarks <-
		list(gp1 = c(
			"HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
			"HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE",
			"HALLMARK_COMPLEMENT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
			"HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME", "HALLMARK_BILE_ACID_METABOLISM",
			"HALLMARK_PROTEIN_SECRETION"
		), gp2 = c(
			"HALLMARK_G2M_CHECKPOINT",
			"HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MITOTIC_SPINDLE",
			"HALLMARK_MYC_TARGETS_V2"
		))

	pull_hallmark_genes <- function(hallmark_gene_set){
		hallmark_gene_sets %>%
			dplyr::filter(gs_name == hallmark_gene_set) %>%
			dplyr::pull(gene_symbol)
	}

	subtype_hallmarks <- unlist(subtype_hallmarks) %>%
		set_names(.)

	hallmark_genes <- map(subtype_hallmarks, pull_hallmark_genes)

	for (i in seq_along(hallmark_genes)){
		hallmark_genes[[i]] <- hallmark_genes[[i]][hallmark_genes[[i]] %in% rownames(seu)]
	}

	seu <- AddModuleScore(
		object = seu,
		features = hallmark_genes,
		ctrl = 5,
		name = 'hallmark'
	)

	names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("hallmark", seq(1, length(hallmark_genes))))] <- names(hallmark_genes)

	return(seu)


}

myseus <- list()
output_plots <- list()
mynbs <- list()

# wu ------------------------------
# wu
study = "wu"
merged_metadata = read_csv(path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/scanpy/merged/metadata.csv"))) %>%
  set_names(c("cell", "sample_id", "merged_leiden")) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, ".h5ad")) %>%
  dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
  identity()

# SRR13884240 numbat ------------------------------
sample_id = "SRR13884240"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)

myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)


output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]]

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)



output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)


output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()

# SRR13884241 numbat ------------------------------
sample_id = "SRR13884241"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)

mynbs[[sample_id]] <- readRDS("~/single_cell_projects/resources/wu_et_al_proj/output/numbat/SRR13884240_numbat.rds")

plot_phylo_heatmap(mynbs[[sample_id]])


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()


# SRR13884242 numbat ------------------------------
sample_id = "SRR13884242"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()

# SRR13884243 numbat in progress ------------------------------
sample_id = "SRR13884243"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)

mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
	dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
	tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

# GT_by_cluster_plot <- ggplot(myseus[[sample_id]]@meta.data, aes(x = `gene_snn_res.0.2`, fill = GT_opt)) +
# 	geom_bar(position="fill") +
# 	coord_flip()
#
# cluster_plot <- DimPlot(myseus[[sample_id]], group.by = "gene_snn_res.0.2")
# GT_plot <- DimPlot(myseus[[sample_id]], group.by = "GT_opt")
#
# GT_by_cluster_plot / (cluster_plot + GT_plot)
#
# DimPlot(myseus[[sample_id]], group.by = "type") +    labs(title = sample_id)
#
# # seu <- seu[,!seu$type %in% c("Red Blood Cells")]
#
# DimPlot(myseus[[sample_id]], group.by = "clone_opt")
#
# DimPlot(myseus[[sample_id]], group.by = "gene_snn_res.0.2")

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.2)

phylo_plot_output <- output_plots[[sample_id]][["numbat_phylo"]]$result

plot_variability_at_SCNA(phylo_plot_output)

# SRR13884244 ------------------------------
sample_id = "SRR13884244"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)

# SRR13884245 ------------------------------
sample_id = "SRR13884245"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)

# SRR13884246 numbat ------------------------------
sample_id = "SRR13884246"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]


myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)

mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()

# SRR13884247 ------------------------------
sample_id = "SRR13884247"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)

# SRR13884248 numbat ------------------------------
sample_id = "SRR13884248"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()

# saveRDS(myseus[[sample_id]], seu_path)

# SRR13884249 numbat ------------------------------
sample_id = "SRR13884249"
study = "wu"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()


# yang------------------------------
# yang
study = "yang"
merged_metadata = read_csv(path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/scanpy/merged/metadata.csv"))) %>%
  set_names(c("cell", "sample_id", "merged_leiden")) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, ".h5ad")) %>%
  dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
  identity()

# SRR14800534 numbat ------------------------------
sample_id = "SRR14800534"
study = "yang"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR14800536 numbat ------------------------------
sample_id = "SRR14800536"
study = "yang"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR14800537 numbat ------------------------------
sample_id = "SRR14800537"
study = "yang"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR14800539 numbat ------------------------------
sample_id = "SRR14800539"
study = "yang"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR14800540 numbat ------------------------------
sample_id = "SRR14800540"
study = "yang"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR14800541 numbat ------------------------------
sample_id = "SRR14800541"
study = "yang"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()


# field ------------------------------
study = "field"
merged_metadata = read_csv(path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/scanpy/merged/metadata.csv"))) %>%
  set_names(c("cell", "sample_id", "merged_leiden")) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, ".h5ad")) %>%
  dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
  identity()

# SRR17960481 numbat ------------------------------
sample_id = "SRR17960481"
study = "field"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

merged_metadata_transfer <-
	merged_metadata %>%
	dplyr::filter(sample_id == {{sample_id}}) %>%
	tibble::column_to_rownames("cell") %>%
	identity()

seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR17960482 numbat ------------------------------
sample_id = "SRR17960482"
study = "field"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR17960483 numbat ------------------------------
sample_id = "SRR17960483"
study = "field"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# SRR17960484 numbat ------------------------------
sample_id = "SRR17960484"
study = "field"
seu_path = path(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/seurat/{sample_id}_infercnv_numbat_seu.rds"))
# read_regress_save(seu_path)
myseus[[sample_id]] <- readRDS(seu_path)

# merged_metadata_transfer <-
# 	merged_metadata %>%
# 	dplyr::filter(sample_id == {{sample_id}}) %>%
# 	tibble::column_to_rownames("cell") %>%
# 	identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)



output_plots[[sample_id]][["dimplot"]] <- DimPlot(myseus[[sample_id]], group.by = c("merged_leiden", "gene_snn_res.0.2"))

output_plots[[sample_id]][["merged_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["sample_marker"]] <- plot_markers(myseus[[sample_id]], metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE)

output_plots[[sample_id]][["merged_marker"]] + output_plots[[sample_id]][["sample_marker"]]

myseus[[sample_id]] <- annotate_seu_with_rb_subtype_gene_expression(myseus[[sample_id]])

# saveRDS(myseus[[sample_id]], seu_path)


mynbs[[sample_id]] <- readRDS(glue("~/single_cell_projects/resources/{study}_et_al_proj/output/numbat/{sample_id}_numbat.rds"))

nb_meta <- mynbs[[sample_id]][["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
  dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
  tibble::column_to_rownames("cell")

myseus[[sample_id]] <- Seurat::AddMetaData(myseus[[sample_id]], nb_meta)

myannot = mynbs[[sample_id]]$clone_post[,c("cell", "GT_opt")]

output_plots[[sample_id]][["numbat_phylo"]] <-  safe_plot_phylo(mynbs[[sample_id]], myseus[[sample_id]], myannot, sample_id, clone_bar = FALSE, p_min = 0.9)

output_plots[[sample_id]][["phylo_probability_plot"]] <- ggplotify::as.ggplot(output_plots[[sample_id]][["numbat_phylo"]][["result"]]) / output_plots[[sample_id]][["numbat_phylo"]][["result"]] %>%
  plot_variability_at_SCNA()



# gather objects ------------------------------

pdf("results/wu/phylo_probability_plots.pdf")
transpose(output_plots)[["phylo_probability_plot"]]
dev.off()
browseURL("results/wu/phylo_probability_plots.pdf")

subtype_hallmarks <-
	list(gp1 = c(
		"HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
		"HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE",
		"HALLMARK_COMPLEMENT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
		"HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME", "HALLMARK_BILE_ACID_METABOLISM",
		"HALLMARK_PROTEIN_SECRETION"
	), gp2 = c(
		"HALLMARK_G2M_CHECKPOINT",
		"HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MITOTIC_SPINDLE",
		"HALLMARK_MYC_TARGETS_V2"
	)) %>%
	unlist()

# plot rb subtype expression (liu et al. 2020)
plot_rb_subtype_expression <- function(seu, seu_name, subtype_hallmarks){

	myfeatures <- c("exprs_gp1", "exprs_gp2", c(subtype_hallmarks))
	featureplots <- FeaturePlot(seu, myfeatures, combine = FALSE)

	featureplots <- map(featureplots, ~(.x + labs(subtitle = seu_name))) %>%
		set_names(myfeatures)

}




pdf("results/test0.pdf")
test0$SRR13884242
dev.off()

browseURL("results/test0.pdf")

subtype_expression_plots <- imap(myseus, plot_rb_subtype_expression, subtype_hallmarks) %>%
	purrr::transpose()

pdf("results/subtype_hallmark_plots.pdf")
subtype_expression_plots
dev.off()

browseURL("results/subtype_hallmark_plots.pdf")

marker_plots <- ls(pattern = ".*[0-9]_marker_plot") %>%
	set_names(.) %>%
	map(get) %>%
	identity()

pdf("results/marker_plots.pdf", width = 10)
marker_plots
dev.off()

browseURL("results/marker_plots.pdf")

dimplots <- ls(pattern = ".*[0-9]_dimplot") %>%
	set_names(.) %>%
	map(get) %>%
	identity()


pdf("results/dimplots.pdf", width = 10)
dimplots
dev.off()

browseURL("results/dimplots.pdf")

# markers of interest ------------------------------



checked_cluster_markers <-
	list(
"0" = c("PCLAF", "TYMS"),
"1" = c("RCVRN"),
"2" = c("TOP2A", "NUSAP1", "RRM2"),
"4" = c("UBE2C", "ARL6IP1", "PTTG1"),
"5" = c("EPB41L4A-AS1", "HSPB1", "DNAJB1"),
"6" = c("GNB3", "PDE6H"),
# "7" = c("B2M", "APOE", "FTL", "VIM", "LGALS1"),
"8" = c("HIST1H4C")
# "9" = c("GNGT1", "ROM1", "GNAT1")
)

plot_markers_by_cell_cycle <- plot_cluster_markers_by_cell_type <- function(seu, checked_cluster_markers){
	cluster_plots <- map(checked_cluster_markers, ~VlnPlot(seu, features = .x, group.by = "Phase"))

	cluster_plots <- map2(cluster_plots, names(checked_cluster_markers), ~(.x + labs(subtitle = .y)))

	return(cluster_plots)

}

cell_cycle_plots <- map(myseus, plot_markers_by_cell_cycle, checked_cluster_markers)

pdf("results/cell_cycle_markers.pdf")
for (i in names(checked_cluster_markers)){
	for(j in cell_cycle_plots){
		print(j[i])
	}
}
dev.off()

browseURL("results/cell_cycle_markers.pdf")

# ------------------------------

checked_marker_vector <- map_chr(checked_cluster_markers, 1)

plot_feature_across_seus <- function(seu, seu_title, myvar){

	varplot <- seu %>%
		FeaturePlot(features = myvar)

	dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2")

	phase_plot <- DimPlot(seu, group.by = "Phase")

	(varplot | (phase_plot / dimplot)) +
		plot_annotation(title = seu_title)

	# mypatch = wrap_plots(varplot, dimplot, phase_plot, nrow = 1) +
	# 	plot_annotation(title = seu_title)
}

common_marker_plots <- imap(myseus, plot_feature_across_seus, checked_marker_vector)

pdf("results/common_marker_plots2.pdf", height = 12, width = 16)
common_marker_plots
dev.off()

browseURL("results/common_marker_plots2.pdf")

compplot_feature_and_clusters <- function(seu, feature){
	fp <- FeaturePlot(seu, feature)

	cp <- DimPlot(seu, group.by = "Phase")

	dp1 <- DimPlot(seu, group.by = c("gene_snn_res.0.15", "merged_leiden"))

	mypatch <- wrap_plots(fp, cp, nrow = 1) / dp1

	return(mypatch)
}

test0 <- map(myseus, compplot_feature_and_clusters, "MT-ND3")

pdf("results/test0.pdf")
test0
dev.off()
browseURL("results/test0.pdf")
