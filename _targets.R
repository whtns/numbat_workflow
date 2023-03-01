## Load your packages, e.g. library(targets).
source("./packages.R")
source("./functions.R")
# debug(make_numbat_heatmaps)
# debug(make_numbat_plot_files)
# debug(diffex_cells)
# debug(filter_numbat_cells)
# debug(diffex_by_cluster)
# debug(enrich_diffex_by_cluster)
# debug(enrich_by_cluster)

tar_option_set(memory = "transient", garbage_collection = TRUE)

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

  # target = function_to_make(arg), ## drake style

  # tar_target(target2, function_to_make2(arg)) ## targets style

  # setup ------------------------------

  tar_target(hallmark_gene_sets, msigdbr(species = "Homo sapiens", category = "H")),
  tar_target(cluster_comparisons, list(
    "SRR13884240" = c(2,3,4),
    "SRR13884241" = c(2,3,4),
    "SRR13884242" = c(2,3,4),
    "SRR13884243" = c(2,3),
    "SRR13884245" = c(2),
    "SRR13884247" = c(2,3),
    "SRR13884249" = c(2),
    "SRR14800534" = c(2),
    "SRR14800537" = c(2,3,4),
    "SRR14800539" = c(2),
    "SRR14800540" = c(2,3),
    "SRR14800541" = c(2,3,4),
    "SRR17960481" = c(2,3),
    "SRR17960482" = c(2),
    "SRR17960483" = c(2,3),
    "SRR17960484" = c(2,3,4)
  )),
  tar_target(wu_seus, unlist(map(list(
    "SRR13884240" = "SRR13884240_infercnv_numbat_seu.rds",
    "SRR13884241" = "SRR13884241_infercnv_numbat_seu.rds",
    "SRR13884242" = "SRR13884242_infercnv_numbat_seu.rds",
    "SRR13884243" = "SRR13884243_seu.rds",
    "SRR13884245" = "SRR13884245_seu.rds",
    "SRR13884247" = "SRR13884247_seu.rds",
    "SRR13884249" = "SRR13884249_infercnv_numbat_seu.rds"
  ), ~ paste0("output/seurat/", .x)))),
  tar_target(wu_nbs, unlist(map(list(
    "SRR13884240" = "SRR13884240_numbat.rds",
    "SRR13884241" = "SRR13884241_numbat.rds",
    "SRR13884242" = "SRR13884242_numbat.rds",
    "SRR13884243" = "SRR13884243_numbat.rds",
    "SRR13884245" = "SRR13884245_numbat.rds",
    "SRR13884247" = "SRR13884247_numbat.rds",
    "SRR13884249" = "SRR13884249_numbat.rds"
  ), ~ paste0("output/numbat/", .x)))),
  tar_target(yang_seus, unlist(map(list(
    "SRR14800534" = "SRR14800534_infercnv_numbat_seu.rds",
    "SRR14800537" = "SRR14800537_infercnv_numbat_seu.rds",
    "SRR14800539" = "SRR14800539_infercnv_numbat_seu.rds",
    "SRR14800540" = "SRR14800540_infercnv_numbat_seu.rds",
    "SRR14800541" = "SRR14800541_infercnv_numbat_seu.rds"
  ), ~ paste0("output/seurat/", .x)))),
  tar_target(yang_nbs, unlist(map(list(
    "SRR14800534" = "SRR14800534_numbat.rds",
    "SRR14800537" = "SRR14800537_numbat.rds",
    "SRR14800539" = "SRR14800539_numbat.rds",
    "SRR14800540" = "SRR14800540_numbat.rds",
    "SRR14800541" = "SRR14800541_numbat.rds"
  ), ~ paste0("output/numbat/", .x)))),
  tar_target(field_seus, unlist(map(list(
    "SRR17960481" = "SRR17960481_infercnv_numbat_seu.rds",
    "SRR17960482" = "SRR17960482_infercnv_numbat_seu.rds",
    "SRR17960483" = "SRR17960483_infercnv_numbat_seu.rds",
    "SRR17960484" = "SRR17960484_infercnv_numbat_seu.rds"
  ), ~ paste0("output/seurat/", .x)))),
  tar_target(field_nbs, unlist(map(list(
    "SRR17960481" = "SRR17960481_numbat.rds",
    "SRR17960482" = "SRR17960482_numbat.rds",
    "SRR17960483" = "SRR17960483_numbat.rds",
    "SRR17960484" = "SRR17960484_numbat.rds"
  ), ~ paste0("output/numbat/", .x)))),
  tar_target(myseus, c(
    wu_seus,
    yang_seus,
    field_seus
  )),
  tar_target(mynbs, c(
    wu_nbs,
    yang_nbs,
    field_nbs
  )),

  tar_target(num_clone_log, c(
    # "SRR13884240",
    "SRR13884241" = 3,
    "SRR13884242" = 4,
    "SRR13884243" = 4,
    "SRR13884245" = 3,
    "SRR13884247" = 4,
    "SRR13884249" = 2,
    "SRR14800534" = 2,
    "SRR14800537" = 4,
    # "SRR14800539",
    "SRR14800540" = 4
    # "SRR14800541",
    # "SRR17960481",
    # "SRR17960482",
    # "SRR17960483",
    # "SRR17960484"
  )),


  # wu------------------------------

  tar_target(wu_metadata, get_merged_metadata("output/scanpy/wu_merged/metadata.csv")),

  ## SRR13884240 numbat ------------------------------

  tar_target(SRR13884240_filter_expressions, list(
    'GT_opt %in% c("2a,10b,8a", "2a,10b", "2a") & p_cnv < 0.8 & seg == "2a"',
    'GT_opt %in% c("") & p_cnv > 0.25 & seg == "2a"'
  )),

  tar_target(SRR13884240_plot_files, make_numbat_plot_files("SRR13884240",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884240_heatmaps, make_numbat_heatmaps("SRR13884240",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = wu_metadata)
             ,
  format = "file"),
  tar_target(SRR13884240_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884240_filter_numbat_probability.pdf",
    "SRR13884240",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884240_filter_expressions
  ),
  format = "file"
  ),
  tar_target(SRR13884240_filtered_cells, filter_numbat_cells(
    "SRR13884240",
    myseus,
    mynbs,
    wu_metadata,
    SRR13884240_filter_expressions
  )),
  tar_target(SRR13884240_diffex, diffex_groups(
    "SRR13884240",
    myseus,
    SRR13884240_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884240"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884240_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884240_filtered_cells, "SRR13884240")
  ),
  tar_target(
    SRR13884240_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884240_filtered_cells, "SRR13884240")
  ),

  tar_target(SRR13884240_diffex_by_cluster, diffex_by_cluster("SRR13884240",
                                                              myseus,
                                                              SRR13884240_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884240"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884240_enrich_by_cluster, enrich_by_cluster("SRR13884240",
                                                              myseus,
                                                              SRR13884240_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884240"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884240_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884240",
    myseus,
    SRR13884240_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884240"]],
    group.by = "clone_opt"
  )),

  tar_target(SRR13884240_numbat_plots, collate_numbat_plots(
    "SRR13884240",
    mynbs
  )),

  tar_target(SRR13884240_enrichment_analysis, enrichment_analysis(SRR13884240_diffex)),


  ## SRR13884241 ------------------------------

  tar_target(SRR13884241_filter_expressions, list(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.25 & seg == "2a"'
  )),

  tar_target(SRR13884241_plot_files, make_numbat_plot_files("SRR13884241",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884241_heatmaps, make_numbat_heatmaps("SRR13884241",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884241_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884241_filter_numbat_probability.pdf",
    "SRR13884241",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884241_filter_expressions
  ),
  format = "file"
  ),
  tar_target(SRR13884241_filtered_cells, filter_numbat_cells(
    "SRR13884241",
    myseus,
    mynbs,
    wu_metadata,
    SRR13884241_filter_expressions
  )),
  tar_target(SRR13884241_diffex, diffex_groups(
    "SRR13884241",
    myseus,
    SRR13884241_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884241"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884241_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884241_filtered_cells, "SRR13884241")
  ),
  tar_target(
    SRR13884241_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884241_filtered_cells, "SRR13884241")
  ),

  tar_target(SRR13884241_diffex_by_cluster, diffex_by_cluster("SRR13884241",
                                                              myseus,
                                                              SRR13884241_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884241"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884241_enrich_by_cluster, enrich_by_cluster("SRR13884241",
                                                              myseus,
                                                              SRR13884241_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884241"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884241_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884241",
    myseus,
    SRR13884241_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884241"]],
    group.by = "clone_opt"
  )),

  tar_target(SRR13884241_numbat_plots, collate_numbat_plots(
    "SRR13884241",
    mynbs
  )),

  tar_target(SRR13884241_enrichment_analysis, enrichment_analysis(SRR13884241_diffex)),






  ## SRR13884242 numbat ------------------------------

  tar_target(SRR13884242_filter_expressions, list(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'
  )),

  tar_target(SRR13884242_plot_files, make_numbat_plot_files("SRR13884242",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884242_heatmaps, make_numbat_heatmaps("SRR13884242",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884242_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884242_filter_numbat_probability.pdf",
    "SRR13884242",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884242_filter_expressions
  ),
  format = "file"
  ),
  tar_target(SRR13884242_filtered_cells, filter_numbat_cells(
    "SRR13884242",
    myseus,
    mynbs,
    wu_metadata,
    SRR13884242_filter_expressions
  )),
  tar_target(SRR13884242_diffex, diffex_groups(
    "SRR13884242",
    myseus,
    SRR13884242_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884242"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884242_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884242_filtered_cells, "SRR13884242")
  ),
  tar_target(
    SRR13884242_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884242_filtered_cells, "SRR13884242")
  ),

  tar_target(SRR13884242_diffex_by_cluster, diffex_by_cluster("SRR13884242",
                                                              myseus,
                                                              SRR13884242_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884242"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884242_enrich_by_cluster, enrich_by_cluster("SRR13884242",
                                                              myseus,
                                                              SRR13884242_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884242"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884242_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884242",
    myseus,
    SRR13884242_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884242"]],
    group.by = "clone_opt"
  )),

  tar_target(SRR13884242_numbat_plots, collate_numbat_plots(
    "SRR13884242",
    mynbs
  )),

  tar_target(SRR13884242_enrichment_analysis, enrichment_analysis(SRR13884242_diffex)),




  ## SRR13884243 numbat ------------------------------

  tar_target(SRR13884243_filter_expressions, list(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "16b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"'
  )),

  tar_target(SRR13884243_plot_files, make_numbat_plot_files("SRR13884243",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884243_heatmaps, make_numbat_heatmaps("SRR13884243",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884243_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884243_filter_numbat_probability.pdf",
    "SRR13884243",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884243_filter_expressions
  ),
  format = "file"
  ),
  tar_target(SRR13884243_filtered_cells, filter_numbat_cells(
    "SRR13884243",
    myseus,
    mynbs,
    wu_metadata,
    SRR13884243_filter_expressions
  )),
  tar_target(SRR13884243_diffex, diffex_groups(
    "SRR13884243",
    myseus,
    SRR13884243_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884243"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884243_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884243_filtered_cells, "SRR13884243")
  ),
  tar_target(
    SRR13884243_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884243_filtered_cells, "SRR13884243")
  ),

  tar_target(SRR13884243_diffex_by_cluster, diffex_by_cluster("SRR13884243",
                                                              myseus,
                                                              SRR13884243_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884243"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884243_enrich_by_cluster, enrich_by_cluster("SRR13884243",
                                                              myseus,
                                                              SRR13884243_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884243"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884243_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884243",
    myseus,
    SRR13884243_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884243"]],
    group.by = "clone_opt"
  )),

  tar_target(SRR13884243_numbat_plots, collate_numbat_plots(
    "SRR13884243",
    mynbs
  )),

  tar_target(SRR13884243_enrichment_analysis, enrichment_analysis(SRR13884243_diffex)),

  ## SRR13884245 numbat ------------------------------

  tar_target(SRR13884245_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "16a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16a"'
  )),

  tar_target(SRR13884245_plot_files, make_numbat_plot_files("SRR13884245",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884245_heatmaps, make_numbat_heatmaps("SRR13884245",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884245_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884245_filter_numbat_probability.pdf",
    "SRR13884245",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884245_filter_expressions
  ),
  format = "file"
  ),
  tar_target(SRR13884245_filtered_cells, filter_numbat_cells(
    "SRR13884245",
    myseus,
    mynbs,
    wu_metadata,
    SRR13884245_filter_expressions
  )),
  tar_target(SRR13884245_diffex, diffex_groups(
    "SRR13884245",
    myseus,
    SRR13884245_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884245"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884245_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884245_filtered_cells, "SRR13884245")
  ),
  tar_target(
    SRR13884245_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884245_filtered_cells, "SRR13884245")
  ),

  tar_target(SRR13884245_diffex_by_cluster, diffex_by_cluster("SRR13884245",
                                                              myseus,
                                                              SRR13884245_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884245"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884245_enrich_by_cluster, enrich_by_cluster("SRR13884245",
                                                              myseus,
                                                              SRR13884245_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR13884245"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR13884245_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884245",
    myseus,
    SRR13884245_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884245"]],
    group.by = "clone_opt"
  )),

  tar_target(SRR13884245_numbat_plots, collate_numbat_plots(
    "SRR13884245",
    mynbs
  )),

  tar_target(SRR13884245_enrichment_analysis, enrichment_analysis(SRR13884245_diffex)),


  ## SRR13884247 numbat ------------------------------

  tar_target(SRR13884247_filter_expressions, list(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "6b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "6b"'
  )),
  tar_target(SRR13884247_plot_files, make_numbat_plot_files("SRR13884247",
    myseus,
    mynbs,
    merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884247_heatmaps, make_numbat_heatmaps("SRR13884247",
    myseus,
    mynbs,
    merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884247_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884247_filter_numbat_probability.pdf",
    "SRR13884247",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884247_filter_expressions
  ),
  format = "file"
  ),
  tar_target(SRR13884247_filtered_cells, filter_numbat_cells(
    "SRR13884247",
    myseus,
    mynbs,
    wu_metadata,
    SRR13884247_filter_expressions
  )),
  tar_target(SRR13884247_diffex, diffex_groups(
    "SRR13884247",
    myseus,
    SRR13884247_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884247"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884247_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884247_filtered_cells, "SRR13884247")
  ),
  tar_target(
    SRR13884247_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884247_filtered_cells, "SRR13884247")
  ),

  tar_target(SRR13884247_diffex_by_cluster, diffex_by_cluster("SRR13884247",
    myseus,
    SRR13884247_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884247"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR13884247_enrich_by_cluster, enrich_by_cluster("SRR13884247",
    myseus,
    SRR13884247_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884247"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR13884247_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884247",
    myseus,
    SRR13884247_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884247"]],
    group.by = "clone_opt"
  )),

  tar_target(SRR13884247_numbat_plots, collate_numbat_plots(
    "SRR13884247",
    mynbs
  )),

  tar_target(SRR13884247_enrichment_analysis, enrichment_analysis(SRR13884247_diffex)),

  ## SRR13884249 numbat ------------------------------

  tar_target(SRR13884249_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'
  )),
  tar_target(SRR13884249_plot_files, make_numbat_plot_files("SRR13884249",
    myseus,
    mynbs,
    merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884249_heatmaps, make_numbat_heatmaps("SRR13884249",
    myseus,
    mynbs,
    merged_metadata = wu_metadata
  ),
  format = "file"),
  tar_target(SRR13884249_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR13884249_filter_numbat_probability.pdf",
    "SRR13884249",
    myseus,
    mynbs,
    merged_metadata = wu_metadata,
    SRR13884249_filter_expressions
  )),
  tar_target(SRR13884249_filtered_cells,
             filter_numbat_cells(
               "SRR13884249",
               myseus,
               mynbs,
               wu_metadata,
               SRR13884249_filter_expressions)),

  # sample_id, myseus, mynbs, merged_metadata, myexpressions

  tar_target(SRR13884249_diffex, diffex_groups(
    "SRR13884249",
    myseus,
    SRR13884249_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884249"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR13884249_diffex_by_cluster, diffex_by_cluster("SRR13884249",
    myseus,
    SRR13884249_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884249"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR13884249_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR13884249",
    myseus,
    SRR13884249_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR13884249"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR13884249_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR13884249_filtered_cells, "SRR13884249")
  ),
  tar_target(
    SRR13884249_pcnv_by_reads,
    plot_pcnv_by_reads(SRR13884249_filtered_cells, "SRR13884249")
  ),
  tar_target(SRR13884249_enrichment_analysis, enrichment_analysis(SRR13884249_diffex)),

  tar_target(SRR13884249_numbat_plots, collate_numbat_plots(
    "SRR13884249",
    mynbs
  )),

  # yang ------------------------------

  tar_target(yang_metadata, get_merged_metadata("output/scanpy/yang_merged/metadata.csv")),

  ## SRR14800534 numbat ------------------------------

  tar_target(SRR14800534_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1d"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1d"'
  )),

  tar_target(SRR14800534_plot_files, make_numbat_plot_files("SRR14800534",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800534_heatmaps, make_numbat_heatmaps("SRR14800534",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800534_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR14800534_filter_numbat_probability.pdf",
    "SRR14800534",
    myseus,
    mynbs,
    merged_metadata = yang_metadata,
    SRR14800534_filter_expressions
  )),
  tar_target(SRR14800534_filtered_cells,
             filter_numbat_cells(
               "SRR14800534",
               myseus,
               mynbs,
               yang_metadata,
               SRR14800534_filter_expressions)),

  # sample_id, myseus, mynbs, merged_metadata, myexpressions

  tar_target(SRR14800534_diffex, diffex_groups(
    "SRR14800534",
    myseus,
    SRR14800534_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800534"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR14800534_diffex_by_cluster, diffex_by_cluster("SRR14800534",
                                                              myseus,
                                                              SRR14800534_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR14800534"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR14800534_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR14800534",
    myseus,
    SRR14800534_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800534"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR14800534_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR14800534_filtered_cells, "SRR14800534")
  ),
  tar_target(
    SRR14800534_pcnv_by_reads,
    plot_pcnv_by_reads(SRR14800534_filtered_cells, "SRR14800534")
  ),
  tar_target(SRR14800534_enrichment_analysis, enrichment_analysis(SRR14800534_diffex)),

  tar_target(SRR14800534_numbat_plots, collate_numbat_plots(
    "SRR14800534",
    mynbs
  )),


  ## SRR14800537 numbat ------------------------------

  tar_target(SRR14800537_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  )),

  tar_target(SRR14800537_plot_files, make_numbat_plot_files("SRR14800537",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800537_heatmaps, make_numbat_heatmaps("SRR14800537",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800537_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR14800537_filter_numbat_probability.pdf",
    "SRR14800537",
    myseus,
    mynbs,
    merged_metadata = yang_metadata,
    SRR14800537_filter_expressions
  )),
  tar_target(SRR14800537_filtered_cells,
             filter_numbat_cells(
               "SRR14800537",
               myseus,
               mynbs,
               yang_metadata,
               SRR14800537_filter_expressions)),

  # sample_id, myseus, mynbs, merged_metadata, myexpressions

  tar_target(SRR14800537_diffex, diffex_groups(
    "SRR14800537",
    myseus,
    SRR14800537_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800537"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR14800537_diffex_by_cluster, diffex_by_cluster("SRR14800537",
                                                              myseus,
                                                              SRR14800537_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR14800537"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR14800537_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR14800537",
    myseus,
    SRR14800537_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800537"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR14800537_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR14800537_filtered_cells, "SRR14800537")
  ),
  tar_target(
    SRR14800537_pcnv_by_reads,
    plot_pcnv_by_reads(SRR14800537_filtered_cells, "SRR14800537")
  ),
  tar_target(SRR14800537_enrichment_analysis, enrichment_analysis(SRR14800537_diffex)),

  tar_target(SRR14800537_numbat_plots, collate_numbat_plots(
    "SRR14800537",
    mynbs
  )),

  ## SRR14800539 numbat ------------------------------

  tar_target(SRR14800539_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'
  )),
  tar_target(SRR14800539_plot_files, make_numbat_plot_files("SRR14800539",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800539_heatmaps, make_numbat_heatmaps("SRR14800539",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800539_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR14800539_filter_numbat_probability.pdf",
    "SRR14800539",
    myseus,
    mynbs,
    merged_metadata = yang_metadata,
    SRR14800539_filter_expressions
  )),
  tar_target(SRR14800539_filtered_cells,
             filter_numbat_cells(
               "SRR14800539",
               myseus,
               mynbs,
               yang_metadata,
               SRR14800539_filter_expressions)),

  # sample_id, myseus, mynbs, merged_metadata, myexpressions

  tar_target(SRR14800539_diffex, diffex_groups(
    "SRR14800539",
    myseus,
    SRR14800539_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800539"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR14800539_diffex_by_cluster, diffex_by_cluster("SRR14800539",
                                                              myseus,
                                                              SRR14800539_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR14800539"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR14800539_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR14800539",
    myseus,
    SRR14800539_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800539"]],
    group.by = "clone_opt")
    ),
  tar_target(
    SRR14800539_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR14800539_filtered_cells, "SRR14800539")
  ),
  tar_target(
    SRR14800539_pcnv_by_reads,
    plot_pcnv_by_reads(SRR14800539_filtered_cells, "SRR14800539")
  ),
  tar_target(SRR14800539_enrichment_analysis, enrichment_analysis(SRR14800539_diffex)),

  tar_target(SRR14800539_numbat_plots, collate_numbat_plots(
    "SRR14800539",
    mynbs
  )),


  ## SRR14800540 numbat ------------------------------

  tar_target(SRR14800540_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  )),

  tar_target(SRR14800540_plot_files, make_numbat_plot_files("SRR14800540",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800540_heatmaps, make_numbat_heatmaps("SRR14800540",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800540_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR14800540_filter_numbat_probability.pdf",
    "SRR14800540",
    myseus,
    mynbs,
    merged_metadata = yang_metadata,
    SRR14800540_filter_expressions
  )),
  tar_target(SRR14800540_filtered_cells,
             filter_numbat_cells(
               "SRR14800540",
               myseus,
               mynbs,
               yang_metadata,
               SRR14800540_filter_expressions)),

  # sample_id, myseus, mynbs, merged_metadata, myexpressions

  tar_target(SRR14800540_diffex, diffex_groups(
    "SRR14800540",
    myseus,
    SRR14800540_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800540"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR14800540_diffex_by_cluster, diffex_by_cluster("SRR14800540",
                                                              myseus,
                                                              SRR14800540_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR14800540"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR14800540_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR14800540",
    myseus,
    SRR14800540_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800540"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR14800540_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR14800540_filtered_cells, "SRR14800540")
  ),
  tar_target(
    SRR14800540_pcnv_by_reads,
    plot_pcnv_by_reads(SRR14800540_filtered_cells, "SRR14800540")
  ),
  tar_target(SRR14800540_enrichment_analysis, enrichment_analysis(SRR14800540_diffex)),

  tar_target(SRR14800540_numbat_plots, collate_numbat_plots(
    "SRR14800540",
    mynbs
  )),


  ## SRR14800541 numbat ------------------------------

  tar_target(SRR14800541_filter_expressions, list(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'
  )),
  tar_target(SRR14800541_plot_files, make_numbat_plot_files("SRR14800541",
                                                            myseus,
                                                            mynbs,
                                                            merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800541_heatmaps, make_numbat_heatmaps("SRR14800541",
                                                        myseus,
                                                        mynbs,
                                                        merged_metadata = yang_metadata
  ),
  format = "file"),
  tar_target(SRR14800541_filtered_plot_files, make_filtered_numbat_plots(
    "results/SRR14800541_filter_numbat_probability.pdf",
    "SRR14800541",
    myseus,
    mynbs,
    merged_metadata = yang_metadata,
    SRR14800541_filter_expressions
  )),
  tar_target(SRR14800541_filtered_cells,
             filter_numbat_cells(
               "SRR14800541",
               myseus,
               mynbs,
               yang_metadata,
               SRR14800541_filter_expressions)),

  # sample_id, myseus, mynbs, merged_metadata, myexpressions

  tar_target(SRR14800541_diffex, diffex_groups(
    "SRR14800541",
    myseus,
    SRR14800541_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800541"]],
    group.by = "clone_opt"
  )),
  tar_target(SRR14800541_diffex_by_cluster, diffex_by_cluster("SRR14800541",
                                                              myseus,
                                                              SRR14800541_filtered_cells,
                                                              ident.1 = 1,
                                                              ident.2 = cluster_comparisons[["SRR14800541"]],
                                                              group.by = "clone_opt"
  )),
  tar_target(SRR14800541_enrich_diffex_by_cluster, enrich_diffex_by_cluster(
    "SRR14800541",
    myseus,
    SRR14800541_filtered_cells,
    ident.1 = 1,
    ident.2 = cluster_comparisons[["SRR14800541"]],
    group.by = "clone_opt"
  )),
  tar_target(
    SRR14800541_pcnv_by_nsnp,
    plot_pcnv_by_nsnp(SRR14800541_filtered_cells, "SRR14800541")
  ),
  tar_target(
    SRR14800541_pcnv_by_reads,
    plot_pcnv_by_reads(SRR14800541_filtered_cells, "SRR14800541")
  ),

  tar_target(SRR14800541_numbat_plots, collate_numbat_plots(
    "SRR14800541",
    mynbs
  )),

  tar_target(SRR14800541_enrichment_analysis, enrichment_analysis(SRR14800541_diffex)),

  tar_target(SRR14800541_infercnv, compare_infercnv(myseus, "SRR14800541")),


  # field ------------------------------

  tar_target(field_metadata, get_merged_metadata("~/single_cell_projects/resources/field_et_al_proj/output/scanpy/merged/metadata.csv")),

  # ## SRR17960481 numbat ------------------------------
  #
  # tar_target(SRR17960481_chr1_filter_expressions, list(
  #   'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "1b"',
  #   'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"'
  # )),
  # tar_target(SRR17960481_chr6_filter_expressions, list(
  #   'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "6b"',
  #   'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "6b"'
  # )),
  # tar_target(SRR17960481_plot_files,
  #            make_numbat_plot_files("SRR17960481",
  #
  #                                   "SRR17960481_infercnv_numbat_seu.rds",
  #                                   "SRR17960481_numbat.rds",
  #                                   merged_metadata = field_metadata
  # )),
  # tar_target(SRR17960481_filtered_plot_files, make_filtered_numbat_plots("results/SRR17960481_filter_numbat_probability.pdf", "SRR17960481",
  #                                                                        "SRR17960481_infercnv_numbat_seu.rds",
  #                                                                        "SRR17960481_numbat.rds",
  #                                                                        merged_metadata = field_metadata,
  #                                                                        SRR17960481_chr1_filter_expressions
  # ),
  # format = "file"
  # ),
  # tar_target(SRR17960481_chr1_filtered_cells, filter_numbat_cells(
  #   "SRR17960481",
  #
  #   "SRR17960481_infercnv_numbat_seu.rds",
  #   "SRR17960481_numbat.rds",
  #   field_metadata,
  #   SRR17960481_chr1_filter_expressions
  # )),
  # tar_target(SRR17960481_diffex, diffex_groups("SRR17960481",
  #                                              myseus,
  #                                              SRR17960481_chr1_filtered_cells,
  #                                              ident.1 = 1,
  #                                              ident.2 = cluster_comparisons[["SRR17960481"]],
  #                                             group.by = "clone_opt"
  # )),
  #
  # tar_target(SRR17960481_pcnv_by_nsnp,
  #            plot_pcnv_by_nsnp(SRR17960481_chr1_filtered_cells, "SRR17960481")
  # ),
  #
  # ## SRR17960482 numbat ------------------------------
  #
  # tar_target(SRR17960482_plot_files, make_numbat_plot_files("SRR17960482",
  #
  #                                                           "SRR17960482_infercnv_numbat_seu.rds",
  #                                                           "SRR17960482_numbat.rds",
  #                                                           field_metadata)),
  #
  # ## SRR17960483 numbat ------------------------------
  #
  # tar_target(SRR17960483_plot_files, make_numbat_plot_files("SRR17960483",
  #
  #                                                           "SRR17960483_infercnv_numbat_seu.rds",
  #                                                           "SRR17960483_numbat.rds",
  #                                                           field_metadata)),
  #
  # ## SRR17960484 numbat ------------------------------
  #
  # tar_target(SRR17960484_plot_files, make_numbat_plot_files("SRR17960484",
  #
  #                                                           "SRR17960484_infercnv_numbat_seu.rds",
  #                                                           "SRR17960484_numbat.rds",
  #                                                           field_metadata), format = "file"),

  # gather objects ------------------------------

  #   ## output plots ------------------------------
  #
      tar_target(output_heatmaps,
                 list(
                   "SRR13884240" = SRR13884240_heatmaps,
                   "SRR13884241" = SRR13884241_heatmaps,
                   "SRR13884242" = SRR13884242_heatmaps,
                   "SRR13884243" = SRR13884243_heatmaps,
                   "SRR13884245" = SRR13884245_heatmaps,
                   "SRR13884247" = SRR13884247_heatmaps,
                   "SRR13884249" = SRR13884249_heatmaps,
                   "SRR14800534" = SRR14800534_heatmaps,
                   "SRR14800537" = SRR14800537_heatmaps,
                   "SRR14800539" = SRR14800539_heatmaps,
                   "SRR14800540" = SRR14800540_heatmaps,
                   "SRR14800541" = SRR14800541_heatmaps
                   # "SRR17960481" = SRR17960481_heatmaps,
                   # "SRR17960482" = SRR17960482_heatmaps,
                   # "SRR17960483" = SRR17960483_heatmaps,
                   # "SRR17960484" = SRR17960484_heatmaps
                   )
      ),
  tar_target(output_plots,
             list(
               "SRR13884240" = SRR13884240_plot_files,
               "SRR13884241" = SRR13884241_plot_files,
               "SRR13884242" = SRR13884242_plot_files,
               "SRR13884243" = SRR13884243_plot_files,
               "SRR13884245" = SRR13884245_plot_files,
               "SRR13884247" = SRR13884247_plot_files,
               "SRR13884249" = SRR13884249_plot_files,
               "SRR14800534" = SRR14800534_plot_files,
               "SRR14800537" = SRR14800537_plot_files,
               "SRR14800539" = SRR14800539_plot_files,
               "SRR14800540" = SRR14800540_plot_files,
               "SRR14800541" = SRR14800541_plot_files
               # "SRR17960481" = SRR17960481_plot_files,
               # "SRR17960482" = SRR17960482_plot_files,
               # "SRR17960483" = SRR17960483_plot_files,
               # "SRR17960484" = SRR17960484_plot_files
             )
  ),
  tar_target(numbat_files,
             list(
               "SRR13884240" = c(SRR13884240_numbat_plots),
               "SRR13884241" = c(SRR13884241_numbat_plots),
               "SRR13884242" = c(SRR13884242_numbat_plots),
               "SRR13884243" = c(SRR13884243_numbat_plots),
               "SRR13884245" = c(SRR13884245_numbat_plots),
               "SRR13884247" = c(SRR13884247_numbat_plots),
               "SRR13884249" = c(SRR13884249_numbat_plots),
               "SRR14800534" = c(SRR14800534_numbat_plots),
               "SRR14800537" = c(SRR14800537_numbat_plots),
               "SRR14800539" = c(SRR14800539_numbat_plots),
               "SRR14800540" = c(SRR14800540_numbat_plots),
               "SRR14800541" = c(SRR14800541_numbat_plots)
               # "SRR17960481" = c(SRR17960481_numbat_plots),
               # "SRR17960482" = c(SRR17960482_numbat_plots),
               # "SRR17960483" = c(SRR17960483_numbat_plots),
               # "SRR17960484" = c(SRR17960484_numbat_plots)
             )
  ),
      tar_target(sample_output_files, output_sample_files(output_plots)),
      tar_target(dimplot_files, compact(map(output_plots, "dimplot"))),
      tar_target(dimplot_pdf, qpdf::pdf_combine(dimplot_files, output = "results/dimplots.pdf"), format = "file"),
      tar_target(numbat_probability_files, compact(map(output_heatmaps, "numbat_probability"))),
      tar_target(numbat_pdf, qpdf::pdf_combine(numbat_probability_files, output = "results/numbat_probability.pdf"), format = "file"),
      tar_target(numbat_phylo_probability_files, compact(map(output_heatmaps, "numbat_phylo_probability"))),
      tar_target(numbat_phylo_pdf, qpdf::pdf_combine(numbat_phylo_probability_files, output = "results/numbat_phylo_probability.pdf"), format = "file"),
      tar_target(combined_marker_files, compact(map(output_plots, "combined_marker"))),
      tar_target(combined_marker_pdf, qpdf::pdf_combine(combined_marker_files, output = "results/combined_marker.pdf"), format = "file"),
      tar_target(clone_distribution_files, compact(map(output_plots, "clone_distribution"))),
      tar_target(clone_distribution_pdf, qpdf::pdf_combine(clone_distribution_files, output = "results/clone_distribution.pdf"), format = "file"),

  tar_target(expression_files, compact(map(numbat_files, "numbat"))),
  tar_target(expression_pdf, numbat_pngs_to_pdf(expression_files, "results/numbat_expression.pdf")),
  tar_target(sridhar_expression_files, compact(map(numbat_files, "numbat_sridhar"))),
  tar_target(sridhar_expression_pdf, numbat_pngs_to_pdf(sridhar_expression_files, "results/numbat_expression_sridhar.pdf")),
  # tar_target(bulk_clone_files, compact(map(numbat_files, "final_clones"))),
  # tar_target(bulk_clones_pdf, numbat_pngs_to_pdf(bulk_clone_files, "results/numbat_bulk_clones.pdf")),


  # #
  #   # ## markers of interest ------------------------------
  #   #
  #   # checked_cluster_markers =
  #   #   list(
  #   #     "0" = c("PCLAF", "TYMS"),
  #   #     "1" = c("RCVRN"),
  #   #     "2" = c("TOP2A", "NUSAP1", "RRM2"),
  #   #     "4" = c("UBE2C", "ARL6IP1", "PTTG1"),
  #   #     "5" = c("EPB41L4A-AS1", "HSPB1", "DNAJB1"),
  #   #     "6" = c("GNB3", "PDE6H"),
  #   #     # "7" = c("B2M", "APOE", "FTL", "VIM", "LGALS1"),
  #   #     "8" = c("HIST1H4C")
  #   #     # "9" = c("GNGT1", "ROM1", "GNAT1")
  #   #   )
  #   #
  #   # cell_cycle_plots = map(myseus, plot_markers_by_cell_cycle, checked_cluster_markers),
  #   #
  #   # pdf("results/cell_cycle_markers.pdf")
  #   # for (i in names(checked_cluster_markers)){
  #   #   for(j in cell_cycle_plots){
  #   #     print(j[i])
  #   #   }
  #   # }
  #   #
  #   # # common marker plots 2
  #   #
  #   # checked_marker_vector = map_chr(checked_cluster_markers, 1),
  #   #
  #   # common_marker_plots = imap(myseus, plot_feature_across_seus, checked_marker_vector)
  #   #
  #   # pdf("results/common_marker_plots2.pdf", height = 12, width = 16)
  #   # common_marker_plots
  #   #
)
