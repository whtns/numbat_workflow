## functions

plot_numbat <- function(nb, myseu, myannot, mytitle, ...) {
  # browser()
  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    identity()

  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")

  # mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  num_cols <-  length(unique(nb$clone_post$clone_opt))
  mypal <- scales::hue_pal()(num_cols) %>%
    set_names(seq(num_cols))

  nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = "GT_opt",
    annot_bar_width = 1,
    ...
  ) +
    labs(title = mytitle)
}

plot_numbat_w_phylo <- function(nb, myseu, myannot, mytitle, ...) {
  browser()
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
    show_phylo = TRUE,
    annot_bar_width = 1,
    ...
  ) +
    labs(title = mytitle)
}

safe_plot_numbat <- safely(plot_numbat, otherwise = NA_real_)

safe_plot_numbat_w_phylo <- safely(plot_numbat_w_phylo, otherwise = NA_real_)

# plot variability at SCNA
plot_variability_at_SCNA <- function(phylo_plot_output, chrom = "1", p_min = 0.9){
  # browser()
  test0 <-
    phylo_plot_output %>%
    dplyr::mutate(seg = factor(seg, levels = str_sort(unique(seg), numeric = TRUE))) %>%
    # dplyr::filter(CHROM == chrom) %>%
    identity()

  p_cnv_plot <- ggplot(test0, aes(x = cell_index, y = p_cnv, color = cnv_state)) +
    geom_point(size = 0.1, alpha = 0.1) +
    geom_hline(aes(yintercept = p_min)) +
    scale_x_reverse() +
    scale_color_manual(values = c("amp" = "#7f180f",
                                  "bamp" = "pink",
                                  "del" = "#010185",
                                  "loh" = "#387229")) +
    facet_wrap(~seg)

  return(p_cnv_plot)
}

read_regress_save <- function(myseus){
  seu <- readRDS(myseus[[sample_id]])

  seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seu))

  # cell cycle effects strongly mitigated in PCA
  seu <- RunPCA(seu, features = VariableFeatures(seu), nfeatures.print = 10)

  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.15)

  seu <- RunUMAP(seu, dims = 1:10, min.dist = 0.01)

  saveRDS(seu, myseus)
}

read_unregress_cc_save <- function(myseus){
  seu <- readRDS(myseus[[sample_id]])

  seu <- seuratTools::clustering_workflow(seu, resolution = c(0.2, 0.4))

  # cell cycle effects strongly mitigated in PCA
  # seu <- seuratTools::seurat_reduce_dimensions(seu)

  # seu <- RunPCA(seu, features = VariableFeatures(seu), nfeatures.print = 10)
  #
  # seu <- RunUMAP(seu, dims = 1:30)

  # saveRDS(seu, myseus)
}

annotate_seu_with_rb_subtype_gene_expression <- function(seu){

  hallmark_gene_sets = msigdbr(species = "Homo sapiens", category = "H")

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

plot_distribution_of_clones_across_clusters <- function(seu, seu_name){
  test0 <- seu@meta.data %>%
    dplyr::mutate(clone_opt = factor(clone_opt))

  clone_plot <- ggplot(test0, aes(x = abbreviation, fill = clone_opt)) +
    geom_bar(position = "fill") +
    coord_flip()

  cluster_plot <- ggplot(test0, aes(x = clone_opt, fill = abbreviation)) +
    geom_bar(position = "fill") +
    coord_flip()

  mypatch <- (clone_plot / cluster_plot) + plot_annotation(title = seu_name)

  return(mypatch)
}

make_numbat_plot_files <- function(done_file, cluster_dictionary){
  # browser()

  output_plots <- list()

  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)


  # merged_metadata_transfer <-
  #   merged_metadata %>%
  #   dplyr::filter(sample_id == {{sample_id}}) %>%
  #   dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
  #   # dplyr::filter(cell %in% colnames(seu)) %>%
  #   tibble::column_to_rownames("cell") %>%
  #   identity()
  #
  # seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)
  #
  # plot_markers(seu, metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE) +
  #   labs(title = sample_id)
  # ggsave(glue("results/{numbat_dir}/{sample_id}_merged_marker.png"))

  plot_markers(seu, metavar = "abbreviation", marker_method = "presto", return_plotly = FALSE) +
    ggplot2::scale_y_discrete(position = "left") +
    labs(title = sample_id)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_sample_marker.pdf"))

  # output_plots[["merged_marker"]] + output_plots[["sample_marker"]]
  # ggsave(glue("results/{numbat_dir}/{sample_id}_combined_marker.pdf"))

  # seu <- annotate_seu_with_rb_subtype_gene_expression(seu)

  DimPlot(seu, group.by = c("abbreviation", "clone_opt", "Phase")) +
    plot_annotation(title = sample_id)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_dimplot.pdf"), width = 10)


  ## clone distribution ------------------------------
  plot_distribution_of_clones_across_clusters(seu, sample_id)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_clone_distribution.pdf"), width = 4, height = 4)

  ## clone tree ------------------------------

  nclones <- length(unique(seu$clone_opt))

  mypal <- scales::hue_pal()(nclones) %>%
    set_names(1:nclones)

  mynb$plot_mut_history(pal = mypal) +
    labs(title = sample_id)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_tree.pdf"), width = 5, height = 2)

  # plot types ------------------------------
  plot_types <- c("dimplot", "sample_marker", "clone_distribution", "tree")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
}

make_numbat_heatmaps <- function(done_file, p_min = 0.9, line_width = 0.1){
  browser()

  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}/"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- mynb$clone_post[,c("cell", "GT_opt")]
  # myannot <- mynb$clone_post[,c("cell", "GT_opt", "clone_opt")]

  ## numbat ------------------------------
  numbat_heatmap <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = TRUE, p_min = p_min, line_width = line_width)[["result"]]
  # numbat_heatmap <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = TRUE, p_min = p_min, line_width = line_width)[["result"]][[1]]

  scna_variability_plot <- plot_variability_at_SCNA(numbat_heatmap[[4]][["data"]])

  patchwork::wrap_plots(numbat_heatmap, scna_variability_plot, ncol = 1)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_probability.pdf"))

  # ## numbat phylo ------------------------------
  # numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)[["result"]]
  #
  # scna_variability_plot_w_phylo <- plot_variability_at_SCNA(numbat_heatmap_w_phylo[[3]][["data"]])
  #
  # patchwork::wrap_plots(numbat_heatmap_w_phylo, scna_variability_plot_w_phylo, ncol = 1)
  # ggsave(glue("results/{numbat_dir}/{sample_id}_numbat_phylo_probability.pdf"))

  plot_types <- c("numbat_probability")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
  # return(numbat_heatmap)
}

make_numbat_heatmaps_old <- function(sample_id, myseus, mynbs, merged_metadata){

  output_plots <- list()
  seu <- readRDS(myseus[[sample_id]])

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{sample_id}}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- mynb$clone_post[,c("cell", "GT_opt", "clone_opt")]

  ## numbat ------------------------------
  output_plots[["numbat"]] <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)[["result"]]

  scna_variability_plot <- plot_variability_at_SCNA(output_plots[["numbat"]][[3]][["data"]])

  patchwork::wrap_plots(output_plots[["numbat"]], scna_variability_plot, ncol = 1)
  ggsave(glue("results/{sample_id}_numbat_probability.pdf"))

  ## numbat phylo ------------------------------
  output_plots[["numbat_phylo"]] <- safe_plot_numbat_w_phylo(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)[["result"]]

  scna_variability_plot <- plot_variability_at_SCNA(output_plots[["numbat_phylo"]][[3]][["data"]])

  patchwork::wrap_plots(output_plots[["numbat_phylo"]], scna_variability_plot, ncol = 1)
  ggsave(glue("results/{sample_id}_numbat_phylo_probability.pdf"))

  plot_types <- c("numbat_phylo_probability","numbat_probability")

  plot_files <- glue("results/{sample_id}_{plot_types}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
}



make_filtered_numbat_plots <- function(output_file, sample_id, myseus, mynbs, merged_metadata, myexpressions){

  output_plots <- list()
  seu <- readRDS(myseus[[sample_id]])

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{sample_id}}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot = mynb$clone_post[,c("cell", "GT_opt", "clone_opt")]

  if(!'' %in% unlist(myexpressions)){
    final_phylo_heatmap <- filter_phylo_plot(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9, expressions = myexpressions)

    test0 <- final_phylo_heatmap[[3]][["data"]] %>%
      plot_variability_at_SCNA()

    patchwork::wrap_plots(final_phylo_heatmap / test0)

    ggsave(output_file)

  }

  return(output_file)
}

filter_numbat_cells <- function(sample_id, myseus, mynbs, merged_metadata, myexpressions){

  seu <- readRDS(myseus[[sample_id]])

  seu_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell")

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{sample_id}}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot = seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, GT_opt, clone_opt, nCount_gene, nFeature_gene) %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-"))

  if(!'' %in% unlist(myexpressions)){
    filtered_nb <- filter_phylo_plot(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9, expressions = myexpressions)
  }

  returned_meta <- dplyr::left_join(
    filtered_nb[[3]][["data"]],
    myannot,
    by = "cell") %>%
    dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
    dplyr::filter(!is.na(cell)) %>%
    identity()

  return(returned_meta)
}

plot_rb_subtype_expression <- function(seu, seu_name, subtype_hallmarks){

  myfeatures <- c("exprs_gp1", "exprs_gp2", c(subtype_hallmarks))
  featureplots <- FeaturePlot(seu, myfeatures, combine = FALSE)

  featureplots <- map(featureplots, ~(.x + labs(subtitle = seu_name))) %>%
    set_names(myfeatures)

}

plot_markers_by_cell_cycle <- plot_cluster_markers_by_cell_type <- function(seu, checked_cluster_markers){
  cluster_plots <- map(checked_cluster_markers, ~VlnPlot(seu, features = .x, group.by = "Phase"))

  cluster_plots <- map2(cluster_plots, names(checked_cluster_markers), ~(.x + labs(subtitle = .y)))

  return(cluster_plots)

}

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

compplot_feature_and_clusters <- function(seu, feature){
  fp <- FeaturePlot(seu, feature)

  cp <- DimPlot(seu, group.by = "Phase")

  dp1 <- DimPlot(seu, group.by = c("gene_snn_res.0.15", "merged_leiden"))

  mypatch <- wrap_plots(fp, cp, nrow = 1) / dp1

  return(mypatch)
}

get_merged_metadata <- function(meta_path){
  metadata = read_csv(meta_path) %>%
    set_names(c("cell", "sample_id", "merged_leiden")) %>%
    dplyr::mutate(sample_id = str_remove(sample_id, ".h5ad")) %>%
    dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
    identity()

  return(metadata)
}

filter_phylo_plot <- function(nb, myseu, myannot, mytitle, expressions, ...) {
  # browser()
  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    identity()

  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")

  mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  initial_phylo_heatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = "GT_opt",
    ...
  ) +
    labs(title = mytitle)

  phylo_heatmap_data <- initial_phylo_heatmap$data %>%
    dplyr::left_join(myannot, by = "cell")

  pull_cells_matching_expression <- function(myexpression, phylo_heatmap_data){
    # browser()
    excluded_cells <-
      phylo_heatmap_data %>%
      dplyr::filter(!!parse_expr(myexpression)) %>%
      dplyr::pull(cell) %>%
      identity()

    return(excluded_cells)

  }

  excluded_cells <- map(expressions, pull_cells_matching_expression, phylo_heatmap_data) %>%
    unlist()

  filtered_myannot <-
    myannot %>%
    dplyr::filter(!cell %in% excluded_cells)

  final_phylo_heatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = filtered_myannot,
    show_phylo = FALSE,
    sort_by = "GT_opt",
    ...
  ) +
    labs(title = mytitle)


  return(final_phylo_heatmap)

}

diffex_groups_old <- function(sample_id, myseus, celldf, ...){

  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  Seurat::FindMarkers(seu, ...)

}


retrieve_numbat_seurat <- function(done_file){

  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  return(seu)

}


diffex_groups <- function(done_file, filter_expressions, idents = NULL){

  sample_id <- path_file(path_dir(done_file))

  ident.1 = idents[[sample_id]]

  filter_expressions <- filter_expressions[[sample_id]]

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  pull_cells_matching_expression <- function(myexpression, joint_post){
    # browser()
    excluded_cells <-
      joint_post %>%
      dplyr::filter(!!parse_expr(myexpression)) %>%
      dplyr::pull(cell) %>%
      identity()

    return(excluded_cells)

  }

  joined_nb_meta <-
    mynb$clone_post %>%
    dplyr::select(-p_cnv) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  excluded_cells <- map(filter_expressions, pull_cells_matching_expression, joined_nb_meta) %>%
    unlist()

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    dplyr::filter(!cell %in% excluded_cells) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  diffex <- Seurat::FindMarkers(seu, ident.1)

  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_gsea.pdf")
  gse_plot <- enrichment_analysis(diffex) +
    labs(title = glue("{sample_id}"))
  ggsave(gse_plot_path)

  diffex_path <- glue("results/{numbat_dir}/{sample_id}_diffex.csv")
  diffex <-
    diffex %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::distinct(ensgene, .keep_all = TRUE)
  write_csv(diffex, diffex_path)

  return(list(sample_id, diffex_path, gse_plot_path))

}

diffex_cells <- function(myseus, cells.1, cells.2, ...){

  seu <- readRDS(myseus[[sample_id]])

  Seurat::FindMarkers(seu$gene, cells.1 = cells.1, cells.2 = cells.2)

}

split_by_cnv <- function(returned_meta, myseg = "2a"){
  test0 <-
    returned_meta %>%
    dplyr::filter(!is.na(cell)) %>%
    dplyr::filter(seg == myseg) %>%
    dplyr::mutate(cnv_status = ifelse(p_cnv > 0.5, "present", "absent")) %>%
    dplyr::group_by(cnv_status) %>%
    dplyr::group_split() %>%
    map(pull, cell) %>%
    identity()
}

plot_pcnv_by_nsnp <- function(tbl, sample_id){
  ggplot(tbl, aes(x = p_cnv, y = n_snp)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_wrap(~seg, ncol = 1)

  ggsave(glue("results/{sample_id}_pcnv_by_nsnp.pdf"))

  return(glue("results/{sample_id}_pcnv_by_nsnp.pdf"))

}

enrichment_analysis <- function(df, fold_change_col = "avg_log2FC", pvalueCutoff = 0.1, ...){
  # browser()

  df <-
    df %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::distinct(entrez, .keep_all = TRUE)

  # we want the log2 fold change
  original_gene_list <- df[[fold_change_col]]

  # name the vector
  names(original_gene_list) <- df$entrez

  # omit any NA values
  gene_list<-na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)

  gse <- clusterProfiler::gseGO(geneList=gene_list,
               ont ="BP",
               keyType = "ENTREZID",
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = pvalueCutoff,
               verbose = TRUE,
               OrgDb = "org.Hs.eg.db",
               pAdjustMethod = "BH")

  gse <- clusterProfiler::simplify(gse)

  possible_dotplot <- possibly(dotplot, otherwise = NA_real_)

  mydotplot <- possible_dotplot(gse, showCategory=10, split=".sign", font.size = 12, ...) + facet_grid(.~.sign)

  if(is.null(mydotplot)){
    mydotplot <- ggplot()
  }

  return(mydotplot)

}

make_cell_cycle_plot <- function(sample_id, myseus){

  output_plots <- list()
  seu <- readRDS(myseus[[sample_id]])

  DimPlot(seu, group.by = c("Phase")) +
    plot_annotation(title = sample_id)
  ggsave(glue("results/{sample_id}_cell_cycle.pdf"), width = 7, height = 7)

  return(glue("results/{sample_id}_cell_cycle.pdf"))
}

diffex_by_cluster <- function(sample_id, myseus, celldf, ...){
  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  clusters <- unique(seu$gene_snn_res.0.2) %>%
    set_names(.)

  filter_seu_to_cluster <- function(cluster, seu){
    seu[,(seu@meta.data$gene_snn_res.0.2 == cluster)]
  }

  split_seu <- map(clusters, filter_seu_to_cluster, seu)

  safe_FindMarkers <- purrr::safely(FindMarkers, otherwise = NA_real_)

  cluster_diffex <- map(split_seu, safe_FindMarkers) %>%
    map("result") %>%
    identity()

  cluster_diffex <- cluster_diffex[!is.na(cluster_diffex)] %>%
    map(tibble::rownames_to_column, "symbol")

  write_xlsx(cluster_diffex, glue("results/{sample_id}_cluster_diffex.xlsx"))

  return(glue("results/{sample_id}_cluster_diffex.xlsx"))

}

enrich_diffex_by_cluster <- function(sample_id, myseus, celldf, ...){
  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  clusters <- unique(seu$gene_snn_res.0.2) %>%
    set_names(.)

  filter_seu_to_cluster <- function(cluster, seu){
    seu[,(seu@meta.data$gene_snn_res.0.2 == cluster)]
  }

  clusters <- janitor::tabyl(seu@meta.data, gene_snn_res.0.2, clone_opt) %>%
    rowwise() %>%
    mutate(empty_clone = min(c_across(!any_of("gene_snn_res.0.2")))) %>%
    dplyr::filter(!empty_clone < 3) %>%
    dplyr::pull(`gene_snn_res.0.2`) %>%
    set_names(.) %>%
    identity()

  split_seu <- map(clusters, filter_seu_to_cluster, seu)

  # cluster_diffex <- map(split_seu, Seurat::FindMarkers, ...)

  possible_FindMarkers <- purrr::possibly(FindMarkers, otherwise = NA_real_)

  cluster_diffex <- map(split_seu, possible_FindMarkers, ...) %>%
    # map("result") %>%
    identity()

  cluster_diffex <- cluster_diffex[!is.na(cluster_diffex)]

  safe_enrichment_analysis <- purrr::safely(enrichment_analysis, otherwise = NA_real_)

  enrich_plots <- map(cluster_diffex, safe_enrichment_analysis) %>%
    map("result") %>%
    identity()

  enrich_plots <- enrich_plots[!is.na(enrich_plots)] %>%
    compact()

  enrich_plots <- compact(enrich_plots) %>%
    imap(~(.x + labs(title = .y))) %>%
    identity()

  pdf(glue("results/{sample_id}_diffex_cluster_enrichment.pdf"), width = 10)
  print(enrich_plots)
  dev.off()

  return(glue("results/{sample_id}_diffex_cluster_enrichment.pdf"))


}

enrich_by_cluster  <- function(done_file){
  # browser()

  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]


  cluster_diffex <- seu@misc$markers$gene_snn_res.0.2$presto %>%
    split(.$Cluster) %>%
    map(tibble::column_to_rownames, "Gene.Name")

  drop_cc_genes <- function(df, cc.genes = Seurat::cc.genes){
    cc_genes <- unlist(cc.genes)

    dplyr::filter(df, !rownames(df) %in% cc_genes)
  }

  cluster_diffex <-
    cluster_diffex %>%
    map(drop_cc_genes)

  safe_enrichment_analysis <- purrr::safely(enrichment_analysis, otherwise = NA_real_)

  enrich_plots <- map(cluster_diffex, safe_enrichment_analysis, fold_change_col = "Average.Log.Fold.Change") %>%
    map("result") %>%
    identity()

  enrich_plots <- compact(enrich_plots) %>%
    imap(~(.x + labs(title = .y))) %>%
    identity()

  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_gsea.pdf")

  pdf(gse_plot_path, width = 8, height = 10)
  print(enrich_plots)
  dev.off()

  return(gse_plot_path)


}

plot_pcnv_by_reads <- function(tbl, sample_id){
  ggplot(tbl, aes(x = p_cnv, y = nCount_gene)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_wrap(~seg, ncol = 1)

  ggsave(glue("results/{sample_id}_pcnv_by_reads.pdf"))

  return(glue("results/{sample_id}_pcnv_by_reads.pdf"))

}

convert_numbat_pngs <- function(done_file){

  numbat_output_dir <- path_dir(done_file)

  sample_id <- fs::path_split(done_file)[[1]][[3]]

  numbat_results_dir <- fs::path_split(done_file)[[1]][[2]]

  numbat_pngs <- dir_ls(numbat_output_dir, glob = "*.png") %>%
    set_names(path_file(.))

  numbat_pdfs <- stringr::str_replace(path_file(numbat_pngs), ".png", ".pdf")

  numbat_pdfs <- glue("results/{numbat_results_dir}/{sample_id}/{numbat_pdfs}")

  numbat_images <- purrr::map(numbat_pngs,image_read) %>%
    imap(~image_annotate(.x, .y, size = 50))

  map2(numbat_images, numbat_pdfs, ~image_write(.x, format = "pdf", .y))

  return(numbat_pdfs)
}

compare_infercnv <- function(myseus, sample_id){
  seu <- readRDS(myseus[[sample_id]])
}

make_all_numbat_plots <- function(numbat_dir, num_iter = 2, min_LLR = 2, genome = "hg38", init_k = 3, gtf = gtf_hg38, overwrite = FALSE){

  sample_id = path_file(numbat_dir)

  print(numbat_dir)

  for (i in seq(num_iter)){

    bulk_clone_path = glue("{numbat_dir}/bulk_clones_{i}.png")
    if(!file.exists(bulk_clone_path) | (overwrite = TRUE)){
      # Plot bulk clones
      bulk_clones <- read_tsv(glue("{numbat_dir}/bulk_clones_{i}.tsv.gz"), col_types = cols())
      p = plot_bulks(bulk_clones, min_LLR = min_LLR, use_pos = TRUE,
                     genome = genome) +
        labs(title = sample_id)

      ggsave(bulk_clone_path, p,
             width = 13, height = 2 * length(unique(bulk_clones$sample)),
             dpi = 250)
      print(glue("plotted {numbat_dir}/bulk_clones_{i}.png"))
    }


    bulk_subtrees_path = glue("{numbat_dir}/bulk_subtrees_{i}.png")
    if(!file.exists(bulk_subtrees_path) | (overwrite = TRUE)){

      # Plot bulk subtrees
      bulk_subtrees <- read_tsv(glue("{numbat_dir}/bulk_subtrees_{i}.tsv.gz"), col_types = cols())
      p = plot_bulks(bulk_subtrees, min_LLR = min_LLR,
                     use_pos = TRUE, genome = genome) +
        labs(title = sample_id)

      ggsave(glue("{numbat_dir}/bulk_subtrees_{i}.png"),
             p, width = 13, height = 2 * length(unique(bulk_subtrees$sample)),
             dpi = 250)
    }

  }

  final_bulk_clones_path = glue("{numbat_dir}/bulk_clones_final.png")
  if(!file.exists(final_bulk_clones_path) | (overwrite = TRUE)){

    bulk_clones <- read_tsv(glue("{numbat_dir}/bulk_clones_final.tsv.gz"), col_types = cols())
    p = plot_bulks(bulk_clones, min_LLR = min_LLR, use_pos = TRUE,
                   genome = genome) +
      labs(title = sample_id)
    ggsave(final_bulk_clones_path, p, width = 13,
           height = 2 * length(unique(bulk_clones$sample)),
           dpi = 250)

  }

  # exp_clust_path = glue("{numbat_dir}/exp_roll_clust.png")
  # if(!file.exists(exp_clust_path)){
  #
  #   # # Plot single-cell smoothed expression magnitude heatmap
  #   gexp_roll_wide <- read_tsv(glue("{numbat_dir}/gexp_roll_wide.tsv.gz"), col_types = cols()) %>%
  #     column_to_rownames("cell")
  #   hc <- readRDS(glue("{numbat_dir}/hc.rds"))
  #   p = plot_exp_roll(gexp_roll_wide = gexp_roll_wide,
  #                     hc = hc, k = init_k, gtf = gtf, n_sample = 10000)
  #   labs(title = sample_id)
  #   ggsave(exp_clust_path, p,
  #          width = 8, height = 4, dpi = 200)
  #
  # }

  # phylo_heatmap_path = glue("{numbat_dir}/phylo_heatmap.png")
#
#   if(!file.exists(phylo_heatmap_path)){
#
#     # # Plot single-cell CNV calls along with the clonal phylogeny
#     nb <- readRDS(glue("{numbat_dir}_numbat.rds"))
#     mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")
#
#     phylo_heatmap <- nb$plot_phylo_heatmap(
#       pal_clone = mypal,
#       show_phylo = TRUE
#     ) +
#       labs(title = sample_id)
#     ggsave(phylo_heatmap_path, width = 13,
#            height = 10,
#            dpi = 250)
#
#   }


    return("success!")

}

retrieve_done_files <- function(numbat_dir, kept_samples = NULL){

  done_files <- fs::dir_ls(numbat_dir, glob = "*done.txt", recurse = TRUE) %>%
    set_names(path_file(path_dir(.)))

  sample_ids <- fs::path_split(done_files) %>%
    map_chr(3)

  if(!is.null(kept_samples)){
    done_files <- done_files[sample_ids %in% kept_samples]
  }

  return(done_files)
}


retrieve_numbat_plot_type <- function(numbat_mini_plots, plot_type = "exp_roll_clust.png"){
  # browser()
  retrieved_plot_types <- map(numbat_mini_plots, ~set_names(.x, fs::path_file(.x))) %>%
    map(plot_type) %>%
    compact() %>%
    unlist() %>%
    set_names(fs::path_file(fs::path_dir(.))) %>%
    identity()


}

reroute_done_to_results_pdf <- function(done_file, label = ""){

  sample_id = path_split(done_file)[[1]][[3]]

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  numbat_pdfs <- dir_ls(glue("results/{numbat_dir}/{sample_id}/"), glob = "*.pdf")

  results_file <- glue("results/{numbat_dir}/{sample_id}{label}.pdf")

  qpdf::pdf_combine(numbat_pdfs, results_file)

  return(results_file)

}

retrieve_snakemake_params <- function(done_file){

  sample_id = path_split(done_file)[[1]][[3]]

  log_file <- path(path_dir(done_file), "log.txt")

  log <- read_lines(log_file)[3:26] %>%
    str_split(" = ") %>%
    transpose() %>%
    identity()

  params <- log[[2]] %>%
    set_names(log[[1]])

  return(list(sample_id, params))

}

retrieve_current_param <- function(current_params, myparam){

  sample_ids <- map(current_params, 1)

  param_values <- map(current_params, 2) %>%
    map(myparam) %>%
    set_names(sample_ids)
}

plot_putative_marker_across_samples <- function(mymarkers, done_files, plot_type = FeaturePlot, group_by = "gene_snn_res.0.2", cluster_dictionary){
  print(mymarkers)
  plot_markers_in_sample <- function(done_file, mymarkers, plot_type = FeaturePlot, group_by = group_by, cluster_dictionary){
    # browser()
    sample_id <- path_file(path_dir(done_file))

    numbat_dir = fs::path_split(done_file)[[1]][[2]]

    dir_create(glue("results/{numbat_dir}"))
    dir_create(glue("results/{numbat_dir}/{sample_id}"))

    seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

    seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

    mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

    nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
      dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
      tibble::column_to_rownames("cell")

    seu <- Seurat::AddMetaData(seu, nb_meta)

    seu <- seu[,!is.na(seu$clone_opt)]

    test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
      dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
      dplyr::select("cell", "abbreviation") %>%
      tibble::column_to_rownames("cell")

    seu <- AddMetaData(seu, test0)

    if(identical(plot_type, VlnPlot)){
      feature_plots <- plot_type(seu, features = mymarkers, group.by = group_by, combine = FALSE, pt.size  = 0.5) %>%
        set_names(mymarkers)
    } else if(identical(plot_type, FeaturePlot)){
      feature_plots <- plot_type(seu, features = mymarkers, combine = FALSE) %>%
        set_names(mymarkers)
    }


    feature_plots <- map(feature_plots, ~(.x + labs(title = sample_id)))

    return(feature_plots)

  }

  sample_ids <- path_file(path_dir(done_files))

  myplots <- map(done_files, plot_markers_in_sample, mymarkers = mymarkers, plot_type = plot_type, group_by = group_by, cluster_dictionary) %>%
    set_names(sample_ids)


  myplots0 <-
    myplots %>%
    transpose() %>%
    imap(~{
      patchwork::wrap_plots(.x) +
        plot_annotation(
          title = .y
        )

    })

  if(identical(plot_type, VlnPlot)){
    plot_type_label = "VlnPlot"
  } else if(identical(plot_type, FeaturePlot)){
    plot_type_label = "FeaturePlot"
  }


  plot_paths <- glue("results/{names(myplots0)}_{plot_type_label}_{group_by}.pdf")

  map2(plot_paths, myplots0, ~ggsave(.x, .y, height = 8, width = 14))

  return(plot_paths)

}

find_clone_diff_in_cluster <- function(done_file){
  # browser()
  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  clone_diff_per_cluster <- function(cluster_for_diffex, seu){
    seu0 <- seu[,seu$gene_snn_res.0.2 == cluster_for_diffex]

    Idents(seu0) <- seu0$clone_opt

    # diffex <- FindAllMarkers(seu0) %>%
    #   dplyr::rename(clone_opt = cluster)

    diffex <- FindMarkers(seu0, group.by = "clone_opt", ident.1 = 1, ident.2 = 2) %>%
      dplyr::mutate(cluster = cluster_for_diffex)

    gse_plot <- enrichment_analysis(diffex) +
      labs(title = glue("{cluster_for_diffex}"))

    return(list("diffex" = diffex, "gse_plot" = gse_plot))

  }

  myclusters <- sort(unique(seu$gene_snn_res.0.2)) %>%
    set_names(.)

  possible_clone_diff_per_cluster <- possibly(clone_diff_per_cluster)

  diffex <- map(myclusters, possible_clone_diff_per_cluster, seu)

  enrich_plots <- map(diffex, 2)

  diffex <- map(diffex, 1) %>%
    compact() %>%
    map(tibble::rownames_to_column, "symbol") %>%
    dplyr::bind_rows(.id = "cluster") %>%
    # dplyr::rename(symbol = gene) %>%
    dplyr::mutate(sample_id = sample_id) %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::select(symbol, description, everything()) %>%
    dplyr::distinct(symbol, cluster, .keep_all = TRUE) %>%
    dplyr::arrange(cluster, p_val_adj) %>%
    identity()


  enrich_plot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_diffex_gsea.pdf")

  pdf(enrich_plot_path, width = 8, height = 10)
  print(enrich_plots)
  dev.off()

  diffex_path <- glue("results/{numbat_dir}/{sample_id}_cluster_diffex.csv")
  write_csv(diffex, diffex_path)

  return(list("enrich_plot" = enrich_plot_path, "diffex" = diffex_path))


}


plot_feature_in_seu <- function(done_file, ...){
  sample_id <- path_file(path_dir(done_file))

  numbat_dir = fs::path_split(done_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(glue("output/{numbat_dir}/{sample_id}_numbat.rds"))

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  Seurat::FeaturePlot(seu, ...)
}


collect_clusters_from_seus <- function(done_files){
  # browser()
  gather_clusters <- function(done_file){
    # browser()
    sample_id <- path_file(path_dir(done_file))

    numbat_dir = fs::path_split(done_file)[[1]][[2]]

    seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))
    # asdf
    clusters <- sort(unique(seu$gene_snn_res.0.2))

    tibble::enframe(clusters)

  }

  sample_ids <- path_file(path_dir(done_files))

  names(done_files) <- sample_ids

  myclusters <- map(done_files, gather_clusters) %>%
    dplyr::bind_rows(.id = "sample_id")

  return(myclusters)

}

read_cluter_dictionary <- function(cluster_dictionary_path){
  cluster_dictionary <- read_csv(cluster_dictionary_path) %>%
    split(.$sample_id)

  return(cluster_dictionary)
}

make_pdf_montages <- function(mini_plot_files, mini_heatmaps){

  mini_plot_files <- map2(mini_plot_files, mini_heatmaps, c)

  sample_ids <-
    mini_plot_files %>%
    map(1) %>%
    map(fs::path_split) %>%
    map(c(1,3))

  names(mini_plot_files) <- sample_ids

  montage_images <- function(plot_files, sample_id){
    browser()
    plot_images <- magick::image_read(plot_files, density = 600)

    my_montage <- magick::image_montage(plot_images, tile = '3', geometry='800x', shadow = FALSE)

    # my_montage <- image_montage(plot_images, geometry = c('x200+10+10', 'x800+10+10', 'x100+10+10'), tile = '3x', shadow = FALSE)

    montage_path <- glue("results/{sample_id}_montage.pdf")

    image_write(my_montage, format = "pdf", montage_path)

    return(montage_path)
  }

  montage_paths <- imap(mini_plot_files, montage_images)

  return(montage_paths)

}

browse_celltype_expression <- function(sridhar_seu, symbol){

  pdf(glue("~/tmp/{symbol}.pdf"))
  FeaturePlot(sridhar_seu, features = c(glue("{symbol}")), split.by = "CellType_predict", combine = FALSE)
  dev.off()
  browseURL(glue("~/tmp/{symbol}.pdf"))

  return(glue("~/tmp/{symbol}.pdf"))

}
