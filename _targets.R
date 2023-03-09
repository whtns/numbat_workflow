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
    "SRR13884240" = c(2,3,4,5,6,7), #2p
    "SRR13884241" = c(2,3,4,5), # 2p
    "SRR13884242" = c(2,3,4), # 16q
    "SRR13884243" = c(2,3), # 16q and 1q
    "SRR13884244" = c(2,3), # nothing; should be 1q
    "SRR13884245" = c(2), # nothing; 19 i guess
    "SRR13884246" = c(2), # nothing; too complicated
    "SRR13884247" = c(2,3,4), # 6p; 1q gain only affects clone 4
    "SRR13884248" = c(2), # 6p
    "SRR13884249" = c(2,3), #1q; 2p only clone 3
    "SRR14800534" = c(3), #1q only clone 3; maybe reconsider
    "SRR14800535" = c(2,3), # 1q; maybe reconsider
    "SRR14800536" = c(2,3,4), #16q
    "SRR14800537" = c(2,3,4,5), #16q
    "SRR14800539" = c(2), # nothing; 18 i guess
    "SRR14800540" = c(2), #16q; just 2
    "SRR14800541" = c(2), #16q is only in clone 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = c(3), #1q/16q only in clone 3
    "SRR17960480" = c(2), #nothing; only 7p in clone 2
    "SRR17960481" = c(2,3,4), #1q/2p/6p
    "SRR17960482" = c(2), # nothing; too complicated
    "SRR17960483" = c(2,3), # nothing; clonal 1q/16q
    "SRR17960484" = c(4) # 16q only in clone 4; 2,3,4 have 1q; interesting sample
  )),

  tar_target(num_clone_log, list(
    "SRR13884240" = 3,
    "SRR13884241" = 3,
    "SRR13884242" = 4,
    "SRR13884243" = 4,
    "SRR13884244" = 3,
    "SRR13884245" = 3,
    "SRR13884246" = 3,
    "SRR13884247" = 4,
    "SRR13884248" = 3,
    "SRR13884249" = 4,
    "SRR14800534" = 2,
    "SRR14800535" = 3,
    "SRR14800536" = 3,
    "SRR14800537" = 4,
    "SRR14800539" = 3,
    "SRR14800540" = 4,
    "SRR14800541" = 3,
    "SRR14800543" = 3,
    "SRR17960480" = 3,
    "SRR17960481" = 3,
    "SRR17960482" = 3,
    "SRR17960483" = 3,
    "SRR17960484" = 3
  )),


  # wu metadata ------------------------------

  tar_target(wu_metadata, get_merged_metadata("output/scanpy/wu_merged/metadata.csv")),

  # yang metadata ------------------------------

  tar_target(yang_metadata, get_merged_metadata("output/scanpy/yang_merged/metadata.csv")),

  # field metadata ------------------------------

  tar_target(field_metadata, get_merged_metadata("output/scanpy/field_merged/metadata.csv")),

  # # collin metadata ------------------------------
  #
  # tar_target(collin_metadata, get_merged_metadata("output/scanpy/collin_merged/metadata.csv")),

  # filter expressions -----------------------------

  tar_target(filter_expressions, list(
    "SRR13884240" = c(
      'GT_opt %in% c("2a,10b,8a", "2a,10b", "2a") & p_cnv < 0.8 & seg == "2a"',
      'GT_opt %in% c("") & p_cnv > 0.25 & seg == "2a"'),
    "SRR13884241" = c(
      'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "2a"',
      'clone_opt %in% c(1) & p_cnv > 0.25 & seg == "2a"'),
    "SRR13884242" = c(
      'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "2a"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'),
    "SRR13884243" = c(
      'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "16b"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"'),
    "SRR13884245" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "16a"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16a"'),
    "SRR13884247" = c(
      'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "6b"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "6b"'),
    "SRR13884249" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'),
    "SRR14800534" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1d"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1d"'),
    "SRR14800537" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'),
    "SRR14800539" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'),
    "SRR14800540" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'),
    "SRR14800541" = c(
      'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
      'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'),
    "SRR17960481" = c(
      'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "1b"',
      'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"')
  )),

  # mini files ------------------------------
  tar_target(
    mini_plot_files,
    make_numbat_plot_files(mini_done_files),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    mini_gseas,
    diffex_groups(mini_done_files, filter_expressions, cluster_comparisons),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    test_heatmap,
    make_numbat_heatmaps("output/numbat_sridhar_mini/SRR14800543/done.txt", p_min = 0.5, line_width = 0.1)
  ),



  tarchetypes::tar_files(mini_done_files, retrieve_done_files("output/numbat_sridhar_mini/")),

  tar_target(
    mini_numbat_plots,
    collate_numbat_plots(mini_done_files),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(mini_numbat_expression, retrieve_numbat_plot_type(mini_numbat_plots, "exp_roll_clust.png")),
  tar_target(mini_numbat_expression_pdf, numbat_pngs_to_pdf(mini_numbat_expression, "results/mini_numbat_expression.pdf")),

  tar_target(mini_numbat_heatmap_files, retrieve_numbat_plot_type(mini_numbat_plots, "panel_2.png")),
  tar_target(mini_numbat_heatmap_pdf, numbat_pngs_to_pdf(mini_numbat_heatmap_files, "results/mini_numbat_heatmaps.pdf")),

  tar_target(mini_numbat_bulk_clones_final, retrieve_numbat_plot_type(mini_numbat_plots, "bulk_clones_final.png")),
  tar_target(mini_numbat_bulk_clones_final_pdf, numbat_pngs_to_pdf(mini_numbat_bulk_clones_final, "results/mini_bulk_clones_final.pdf")),

  tar_target(
    mini_numbat_sample_pdfs,
    reroute_done_to_results_pdf(mini_done_files, "mini"),
    pattern = map(mini_done_files),
    iteration = "list",
    format = "file"
  ),


  # large files ------------------------------

  # tar_target(
  #   large_plot_files,
  #   make_numbat_plot_files(large_done_files),
  #   pattern = map(large_done_files),
  #   iteration = "list"
  # ),

  tarchetypes::tar_files(large_done_files, retrieve_done_files("output/numbat_sridhar/")),

  tar_target(
    large_numbat_plots,
    collate_numbat_plots(large_done_files),
    pattern = map(large_done_files),
    iteration = "list"
  ),

  tar_target(large_numbat_expression, retrieve_numbat_plot_type(large_numbat_plots, "exp_roll_clust.png")),
  tar_target(large_numbat_expression_pdf, numbat_pngs_to_pdf(large_numbat_expression, "results/large_numbat_expression.pdf")),

  tar_target(large_numbat_bulk_clones_final, retrieve_numbat_plot_type(large_numbat_plots, "bulk_clones_final.png")),
  tar_target(large_numbat_bulk_clones_final_pdf, numbat_pngs_to_pdf(large_numbat_bulk_clones_final, "results/large_bulk_clones_final.pdf")),

  tar_target(
    large_numbat_sample_pdfs,
    reroute_done_to_results_pdf(large_done_files, "large"),
    pattern = map(large_done_files),
    iteration = "list"
  ),


)

# notes 2023-03-07 ------------------------------

  # yes ------------------------------

  # SRR17960481: 6p interesting; cluster 2 and 4; cluster 2 has mito genes

  # SRR17960483: cluster 6 maybe interesting

  # SRR17960484: cluster 1 enriched for wt clone 1; also has high Xist expression

  # SRR14800541: clear 6p and 16q clones

  # SRR14800540: clear 16q clone

  # SRR14800535: 16q clone

  # SRR13884247: likely 16q clone; no GSEA diffex output

  # SRR13884242: likely 16q clone

  # maybe ------------------------------
  # SRR17960480: possible minor 16q- clone

  # SRR14800543: possible minor 16q/1q clone

  # SRR14800539: possible 16q clone

  # SRR14800537: possible 16q clone

  # SRR14800536: possible 16q clone

  # SRR14800534: possible 16q clone

  # SRR13884249: possible 1q/2p/16q clone

  # SRR13884248: possible 2p

  # SRR13884246: possible 16q clone

  # SRR13884245: possible 1q clone

  # SRR13884244: possible 1q clone

  # SRR13884243: possible 1q clone

  # SRR13884241: possible 1q clone

  # SRR13884240: possible 2p clone

  # no ------------------------------
  # SRR17960482: too complicated



