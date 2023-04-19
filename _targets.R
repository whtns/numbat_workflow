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
# debug(find_diffex_bw_clusters_for_each_clone)

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

  tar_target(cluster_dictionary, read_cluster_dictionary("data/cluster_dictionary.csv")),

  tar_target(hallmark_gene_sets, msigdbr(species = "Homo sapiens", category = "H")),
  tar_target(cluster_comparisons, list(
    "SRR13884240" = c(2,3,4,5,6,7), #2p
    "SRR13884241" = c(2,3,4,5), # 2p
    "SRR13884242" = c(2,3,4), # 16q
    "SRR13884243" = c(2,3), # 16q and 1q
    "SRR13884244" = c(2,3), # nothing; should be 1q
    "SRR13884245" = c(2), # nothing; 19 i guess
    "SRR13884246" = c(2), # nothing; too complicated
    "SRR13884247" = c(2,3,4), # 6p; 1q gain only affects GT 4
    "SRR13884248" = c(2), # 6p
    "SRR13884249" = c(2,3), #1q; 2p only GT 3
    "SRR14800534" = c(3), #1q only GT 3; maybe reconsider
    "SRR14800535" = c(2,3), # 1q; maybe reconsider
    "SRR14800536" = c(2,3,4), #16q
    "SRR14800537" = c(2,3,4,5), #16q
    "SRR14800539" = c(2), # nothing; 18 i guess
    "SRR14800540" = c(2), #16q; just 2
    "SRR14800541" = c(2), #16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = c(3), #1q/16q only in GT 3
    "SRR17960480" = c(2), #nothing; only 7p in GT 2
    "SRR17960481" = c(2,3,4), #1q/2p/6p
    "SRR17960482" = c(2), # nothing; too complicated
    "SRR17960483" = c(2,3), # nothing; clonal 1q/16q
    "SRR17960484" = c(4) # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),

  tar_target(large_clone_comparisons, list(
    "SRR13884242" = list("1_v_2_1q+_16q-" = c("1c", "16b"), "2_v_3_8p+_11p+" = c("8a", "11b"), "3_v_4_2+" = c("2c")), # 16q
    "SRR13884243" = list("1_v_2_1q+_16q-"= c("1b", "16b"), "2_v_3_8p+_11p+" = c("8a", "11a")), # 16q and 1q
    "SRR13884247" = list("1_v_2_6p+" = c("6c"), "2_v_3_17q+" = c("17d"), "3_v_4_10q+" = c("10b"), "3_v_5_1p-" = c("1a"), "2_v_6_2p+" = c("2b")), # 6p; 1q gain only affects GT 4
    "SRR13884249" = list("1_v_2_1q+" = c("1b"), "2_v_3_2p+" = c("2a")), #1q; 2p only GT 3
    "SRR14800534" = list("1_v_2_16q-" = c("16c"), "2_v_3_1q+" = c("1b")), #1q only GT 3; maybe reconsider
    "SRR14800535" = list("1_v_2_16q-" = c("16b"), "2_v_3_1q+" = c("1b")), #1q only GT 3; maybe reconsider
    "SRR14800536" = list("1_v_2_16q-" = c("16b"), "2_v_3_1q+" = c("1b"), "3_v_4_15q+" = c("15b")), #16q
    "SRR14800540" = list("1_v_2_6p+_16q-" = c("6b", "16e"), "2_v_3_1q+_7q+" = c("1d", "1e", "7b")), #16q; just 2
    "SRR14800541" = list("1_v_2_1q+_6p+" = c("1f", "6a"), "2_v_3_16q-" = c("16d")), #16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = list("1_v_2_13loh" = c("13a"), "2_v_3_1q+_16q-" = c("1c", "16c")), #1q/16q only in GT 3
    "SRR17960481" = list("1_v_2_1q+_6p+" = c("1b", "6a", "6b"), "2_v_3_2p+" = c("2a")), #1q/2p/6p
    "SRR17960484" = list("1_v_2_1q+" = c("1c"), "2_v_3_6p+" = c("6a"), "3_v_4_2p+" = c("2a")) # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),

  tar_target(mini_clone_comparisons, list(
    "SRR13884242" = list("1_v_2+_1q+_16q-" = c("1b", "16b"), "2_v_3_8p+_11p+" = c("8a", "11a")), # 16q
    "SRR13884243" = list("1_v_2+_1q+_16q-"= c("1b", "16b"), "2_v_3_8p+_11p+" = c("8a", "11a")), # 16q and 1q
    "SRR13884247" = list("1_v_2+_6p+" = c("6c"), "2_v_3_1q+" = c("1b"), "3_v_4_20p-" = c("20d")), # 6p; 1q gain only affects GT 4
    "SRR13884249" = list("1_v_2+_1q+" = c("1b"), "2_v_3_2p+" = c("2a")), #1q; 2p only GT 3
    "SRR14800534" = list("1_v_2+_16q-" = c("16c"), "2_v_3_1q+" = c("1b")), #1q only GT 3; maybe reconsider
    "SRR14800535" = list("1_v_2+_16q-" = c("16b")), #1q only GT 3; maybe reconsider
    "SRR14800536" = list("1_v_2+_1q+" = c("1b"), "2_v_3_11q-" = c("11b"), "2_v_4_15q+, 19q-" = c("19b", "15b")), #16q
    "SRR14800540" = list("1_v_2+_1q+" = c("1e"), "1_v_5_6p+" = c("6b"), "2_v_3_6q+_13p-_16q-_20q+" = c("6d", "13d", "16b", "16f", "20c"), "2_v_4_3q-" = c("3b")), #16q; just 2
    "SRR14800541" = list("1_v_2_1q+_6p+" = c("1f", "6a"), "2_v_3_16q-" = c("16b"), "2_v_4_8p-" = c("8a"), "4_v_5_19p+" = c("19a")), #16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = list("1_v_2_13loh_18+" = c("13b", "18a"), "2_v_3_1q+" = c("1c"), "3_v_4_16q-_17p-" = c("16c", "17a")), #1q/16q only in GT 3
    "SRR17960481" = list("1_v_2_1q+, 6p+" = c("1b", "6a", "6b"), "2_v_3_2p+_9q+" = c("2a", "9b")), #1q/2p/6p
    "SRR17960484" = list("1_v_2_1q+, 6p+" = c("1e", "6a"), "2_v_3_2p+_10q+" = c("2a", "10b"), "3_v_4_15+" = c("15b"), "2_v_5_16q-" = c("16b"), "5_v_6_17q+" = c("17b"), "5_v_7_5q-" = c("5d"))
    )
    ),# 16q only in GT 4; 2,3,4 have 1q; interesting sample)),

  tar_target(cluster_for_diffex, list(
    "SRR13884242" = 3, # 16q
    "SRR13884243" = 4, # 16q and 1q
    "SRR13884247" = 5, # 6p; 1q gain only affects GT 4
    "SRR13884249" = 2, #1q; 2p only GT 3
    "SRR14800534" = 2, #1q only GT 3; maybe reconsider
    "SRR14800535" = 1, # 1q; maybe reconsider
    "SRR14800536" = 2, #16q
    "SRR14800540" = 4, #16q; just 2
    "SRR14800541" = 1, #16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = 2, #1q/16q only in GT 3
    "SRR17960481" = 2, #1q/2p/6p
    "SRR17960484" = 1 # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),

  tar_target(future_k_init, list(
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

  tar_target(
    current_params,
    retrieve_snakemake_params(mini_done_files),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    current_init_k,
    retrieve_current_param(current_params, "init_k")
  ),

  # wu metadata ------------------------------

  tar_target(wu_metadata, get_merged_metadata("output/scanpy/wu_merged/metadata.csv")),

  # yang metadata ------------------------------

  tar_target(yang_metadata, get_merged_metadata("output/scanpy/yang_merged/metadata.csv")),

  # field metadata ------------------------------

  tar_target(field_metadata, get_merged_metadata("output/scanpy/field_merged/metadata.csv")),

  # # collin metadata ------------------------------
  #
  # tar_target(collin_metadata, get_merged_metadata("output/scanpy/collin_merged/metadata.csv")),


# mini files ------------------------------

tar_target(mini_filter_expressions, list(
  "SRR13884241" = c(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.25 & seg == "2a"'
  ),
  "SRR13884242" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"',
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "16b"'
  ),
  "SRR13884243" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "16b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"',
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "1b"'
  ),
  "SRR13884245" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "16a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16a"'
  ),
  "SRR13884247" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "6c"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "6c"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"',
    'clone_opt %in% c(3,4) & p_cnv < 0.8 & seg == "1b"'
  ),
  "SRR13884249" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "2a"',
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "2a"'
  ),
  "SRR14800534" = c(
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR14800535" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "16b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"'
  ),
  "SRR14800536" = c(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR14800537" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR14800539" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'
  ),
  "SRR14800540" = c(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "1e"',
    'clone_opt %in% c(1,5) & p_cnv > 0.1 & seg == "1e"',
    'clone_opt %in% c(1,2,3) & p_cnv > 0.1 & seg == "6b"',
    'clone_opt %in% c(4,5) & p_cnv < 0.8 & seg == "6b"'
  ),
  "SRR14800541" = c(
    'clone_opt %in% c(2,3,4,5) & p_cnv < 0.8 & seg == "1f"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1f"',
    'clone_opt %in% c(1,2,4,5) & p_cnv > 0.1 & seg == "16b"'
  ),
  "SRR14800543" = c(
    'clone_opt %in% c(3,4) & p_cnv < 0.8 & seg == "1c"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1c"',
    'clone_opt %in% c(1,2,3) & p_cnv > 0.1 & seg == "16c"'
  ),
  "SRR17960481" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR17960484" = c(
    'clone_opt %in% c(2,3,4,5,6,7) & p_cnv < 0.8 & seg == "1e"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1e"',
    'clone_opt %in% c(3,4) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1,2,5,6,7) & p_cnv > 0.1 & seg == "2a"'
    # 'clone_opt %in% c(2,3,4,7) & p_cnv < 0.8 & seg == "6a"',
    # 'clone_opt %in% c(1,5,6) & p_cnv > 0.1 & seg == "6a"',
    # 'clone_opt %in% c(1,2,3,4) & p_cnv > 0.1 & seg == "16b"',
    # 'clone_opt %in% c(5,6,7) & p_cnv < 0.8 & seg == "16b"'
  )
)),

  tar_target(interesting_samples,
             c(
               "SRR13884242",
               "SRR13884243",
               "SRR13884247",
               "SRR13884249",
               "SRR14800534",
               "SRR14800535",
               "SRR14800536",
               "SRR14800540",
               "SRR14800541",
               "SRR14800543",
               "SRR17960481",
               "SRR17960484"
             )),


  tarchetypes::tar_files(mini_done_files, retrieve_done_files("output/numbat_sridhar_mini/", interesting_samples), format = "file"),

  tar_target(
    mini_numbat_pdfs,
    convert_numbat_pngs(mini_done_files),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    mini_plot_files,
    make_numbat_plot_files(mini_done_files, cluster_dictionary),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    mini_dists,
    score_filtration(mini_done_files, cluster_dictionary, mini_filter_expressions),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    filtered_mini_plot_files,
    make_numbat_plot_files(mini_done_files, cluster_dictionary, mini_filter_expressions, extension = "_filtered"),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  # tar_target(
  #   filtered_mini_plot_files,
  #   make_numbat_plot_files("output/numbat_sridhar_mini/SRR14800534/done.txt", cluster_dictionary, mini_filter_expressions),
  # ),

  tar_target(
    mini_heatmaps,
    make_numbat_heatmaps(mini_done_files, p_min = 0.5, line_width = 0.1, mini_filter_expressions, extension = "_filtered"),
    pattern = map(mini_done_files),
    iteration = "list"
  ),


  tar_target(mini_montage_pdfs0, make_pdf_montages(filtered_mini_plot_files, mini_heatmaps)),


  # tar_target(
  #   mini_gseas,
  #   diffex_groups(mini_done_files, mini_filter_expressions, cluster_comparisons),
  #   pattern = map(mini_done_files),
  #   iteration = "list"
  # ),

  # tar_target(
  #   mini_cluster_gseas,
  #   enrich_by_cluster(mini_done_files),
  #   pattern = map(mini_done_files),
  #   iteration = "list"
  # ),

  # in segment ------------------------------
  # tar_target(mini_in_segment_diffex_clones,
  #            find_diffex_clones(mini_done_files, mini_clone_comparisons, location = "in_segment"),
  #            pattern = map(mini_done_files),
  #            iteration = "list"
  # ),
  #
  # tar_target(
  #   mini_in_segment_diffex_clones_for_each_cluster,
  #   find_diffex_bw_clones_for_each_cluster(mini_done_files, mini_clone_comparisons, cluster_dictionary, location = "in_segment"),
  #   pattern = map(mini_done_files),
  #   iteration = "list"
  # ),
  #
  # tar_target(table_mini_in_segment_diffex_clones,
  #            tabulate_diffex_clones(mini_in_segment_diffex_clones_for_each_cluster,
  #                                   "results/diffex_bw_clones_per_cluster_mini_in_segment.xlsx",
  #                                   "results/diffex_bw_clones_per_cluster_mini_in_segment_by_chr.xlsx",
  #                                   mini_in_segment_diffex_clones,
  #                                   "results/diffex_bw_clones_mini_in_segment.xlsx",
  #                                   "results/diffex_bw_clones_mini_in_segment_by_chr.xlsx"),
  # ),

  # out of segment ------------------------------
  # tar_target(mini_out_of_segment_diffex_clones,
  #            find_diffex_clones(mini_done_files, mini_clone_comparisons, location = "out_of_segment"),
  #            pattern = map(mini_done_files),
  #            iteration = "list"
  # ),
  #
  # tar_target(
  #   mini_out_of_segment_diffex_clones_for_each_cluster,
  #   find_diffex_bw_clones_for_each_cluster(mini_done_files, mini_clone_comparisons, cluster_dictionary, location = "out_of_segment"),
  #   pattern = map(mini_done_files),
  #   iteration = "list"
  # ),
  #
  # tar_target(table_mini_out_of_segment_diffex_clones,
  #            tabulate_diffex_clones(mini_out_of_segment_diffex_clones_for_each_cluster,
  #                                   "results/diffex_bw_clones_per_cluster_mini_out_of_segment.xlsx",
  #                                   "results/diffex_bw_clones_per_cluster_mini_out_of_segment_by_chr.xlsx",
  #                                   mini_out_of_segment_diffex_clones,
  #                                   "results/diffex_bw_clones_mini_out_of_segment.xlsx",
  #                                   "results/diffex_bw_clones_mini_out_of_segment_by_chr.xlsx"),
  # ),

  tar_target(
    mini_diffex_bw_clusters_for_each_clone,
    find_diffex_bw_clusters_for_each_clone(mini_done_files, cluster_dictionary),
    pattern = map(mini_done_files),
    iteration = "list"
  ),

  tar_target(
    collated_mini_heatmaps,
    qpdf::pdf_combine(unlist(mini_heatmaps), "results/numbat_sridhar_mini/mini_heatmaps.pdf"),
  ),

  tar_target(mini_numbat_expression, retrieve_numbat_plot_type(mini_numbat_pdfs, "exp_roll_clust.pdf")),
  tar_target(mini_numbat_expression_pdf, qpdf::pdf_combine(mini_numbat_expression, "results/numbat_sridhar_mini/mini_numbat_expression.pdf")),

  tar_target(mini_numbat_heatmap_files, retrieve_numbat_plot_type(mini_numbat_pdfs, "panel_2.pdf")),
  tar_target(mini_numbat_heatmap_pdf, qpdf::pdf_combine(mini_numbat_heatmap_files, "results/numbat_sridhar_mini/mini_numbat_heatmaps.pdf")),

  tar_target(mini_numbat_bulk_clones_final, retrieve_numbat_plot_type(mini_numbat_pdfs, "bulk_clones_final.pdf")),
  tar_target(mini_numbat_bulk_clones_final_pdf, qpdf::pdf_combine(mini_numbat_bulk_clones_final, "results/numbat_sridhar_mini/mini_bulk_clones_final.pdf")),

  # mini_numbat_sample_pdfs ------------------------------
  tar_target(
    mini_numbat_sample_pdfs,
    reroute_done_to_results_pdf(mini_done_files, "_mini"),
    pattern = map(mini_done_files),
    iteration = "list",
    format = "file"
  ),

  tar_target(
    interesting_genes,
    unlist(list(
      "1q" = c("ANP32E", "ASH1L", "ASPM", "CENPF", "CKS1B", "CNIH4", "CRABP2", "ENAH", "ENO1", "IPO9", "KIF14", "NEK2", "NENF", "NUF2", "PSMD4", "RXRG", "UBE2T"),
      "2p" = c("CEP68", "MEIS1", "MDH1", "OST4", "PDIA6", "POLE4", "RAB1A", "SNRPG", "SOX11"),
      "6p" = c("ARID1B", "CASP8AP2", "CLIC1", "CUTA", "DDX39B", "DEK", "DST", "GLO1", "HDAC2", "HDDC2", "HMGA1", "LIN28B", "MCM3", "MNDA", "PRDM1", "SOX4"),
      "16q" = c("CHD9", "CNGB1", "CYBA", "MT1E", "MT1X", "MT2A")
    ))
  ),

  tar_target(
    cluster_markers,
    collect_all_markers(mini_done_files, "results/clusters.xlsx")
  ),

  tar_target(
    marker_gene_featureplots,
    plot_putative_marker_across_samples(interesting_genes, large_done_files, plot_type = FeaturePlot, group_by = "clone_opt", cluster_dictionary),
  ),

  # tar_target(
  #   marker_gene_vlnplots_by_cluster,
  #   plot_putative_marker_across_samples(interesting_genes, mini_done_files, VlnPlot, group_by = "abbreviation", cluster_dictionary),
  # ),
  #
  tar_target(
    marker_gene_vlnplots_by_clone,
    plot_putative_marker_across_samples(interesting_genes, large_done_files, plot_type = VlnPlot, group_by = "clone_opt", cluster_dictionary),
  ),
#
  tar_target(mini_clusters, collect_clusters_from_seus(large_done_files)),

# # large files ------------------------------

tar_target(large_filter_expressions, list(
  "SRR13884241" = c(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.25 & seg == "2a"'
  ),
  "SRR13884242" = c(
    'clone_opt %in% c(2,3,4) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"',
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "16b"'
  ),
  "SRR13884243" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "16b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"',
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "1b"'
  ),
  "SRR13884245" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "16a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16a"'
  ),
  "SRR13884247" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "6c"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "6c"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"',
    'clone_opt %in% c(3,4) & p_cnv < 0.8 & seg == "1b"'
  ),
  "SRR13884249" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "2a"',
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "2a"'
  ),
  "SRR14800534" = c(
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR14800535" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "16b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "16b"'
  ),
  "SRR14800536" = c(
    'clone_opt %in% c(3,4) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR14800537" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR14800539" = c(
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "2a"'
  ),
  "SRR14800540" = c(
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "1e"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1e"',
    'clone_opt %in% c(1,3) & p_cnv > 0.1 & seg == "6b"',
    'clone_opt %in% c(2) & p_cnv < 0.8 & seg == "6b"'
  ),
  "SRR14800541" = c(
    'clone_opt %in% c(2,3,4,5,6,7) & p_cnv < 0.8 & seg == "1f"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1f"',
    'clone_opt %in% c(3,4) & p_cnv < 0.8 & seg == "16b"'
  ),
  "SRR14800543" = c(
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "1c"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "1c"',
    'clone_opt %in% c(1,2,4) & p_cnv > 0.1 & seg == "16c"',
    'clone_opt %in% c(3) & p_cnv < 0.8 & seg == "16c"'

  ),
  "SRR17960481" = c(
    'clone_opt %in% c(2,3) & p_cnv < 0.8 & seg == "1b"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1b"'
  ),
  "SRR17960484" = c(
    'clone_opt %in% c(2,3,4,5) & p_cnv < 0.8 & seg == "1c"',
    'clone_opt %in% c(1) & p_cnv > 0.1 & seg == "1c"',
    'clone_opt %in% c(4,5) & p_cnv < 0.8 & seg == "2a"',
    'clone_opt %in% c(1,2,3) & p_cnv > 0.1 & seg == "2a"',
    'clone_opt %in% c(1,2) & p_cnv > 0.1 & seg == "6a"',
    'clone_opt %in% c(3,4,5) & p_cnv < 0.8 & seg == "6a"'
    # 'clone_opt %in% c(1,5,6) & p_cnv > 0.1 & seg == "6a"',
    # 'clone_opt %in% c(1,2,3,4) & p_cnv > 0.1 & seg == "16b"',
    # 'clone_opt %in% c(5,6,7) & p_cnv < 0.8 & seg == "16b"'
  )
)),

tarchetypes::tar_files(large_done_files, retrieve_done_files("output/numbat_sridhar/", interesting_samples), format = "file"),

tar_target(
  large_numbat_pdfs,
  convert_numbat_pngs(large_done_files),
  pattern = map(large_done_files),
  iteration = "list"
),

tar_target(
  large_plot_files,
  make_numbat_plot_files(large_done_files, cluster_dictionary),
  pattern = map(large_done_files),
  iteration = "list"
),

tar_target(
  filtered_large_plot_files,
  make_numbat_plot_files(large_done_files, cluster_dictionary, large_filter_expressions, extension = "_filtered"),
  pattern = map(large_done_files),
  iteration = "list"
),

tar_target(
  large_heatmaps,
  make_numbat_heatmaps(large_done_files, p_min = 0.5, line_width = 0.1, large_filter_expressions, extension = "_filtered"),
  pattern = map(large_done_files),
  iteration = "list"
),

tar_target(large_montage_pdfs, make_pdf_montages(filtered_large_plot_files, large_heatmaps)),

tar_target(expr_heatmap_pdfs, make_expression_heatmap_comparison(large_numbat_pdfs, large_heatmaps)),

# # tar_target(
# #   large_gseas,
# #   diffex_groups(large_done_files, large_filter_expressions, cluster_comparisons),
# #   pattern = map(large_done_files),
# #   iteration = "list"
# # ),
#
# # tar_target(
# #   large_cluster_gseas,
# #   enrich_by_cluster(large_done_files),
# #   pattern = map(large_done_files),
# #   iteration = "list"
# # ),

  # in segment ------------------------------
  tar_target(large_in_segment_diffex_clones,
             find_diffex_clones(large_done_files, large_clone_comparisons, location = "in_segment"),
             pattern = map(large_done_files),
             iteration = "list"
  ),

  tar_target(
    large_in_segment_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(large_done_files, large_clone_comparisons, cluster_dictionary, location = "in_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),

  tar_target(volcano_thresholds_in_segment,
             list(
               "SRR13884242" = 4,
               "SRR13884243" = 4,
               "SRR13884247" = 4,
               "SRR13884249" = 4,
               "SRR14800534" = 2,
               "SRR14800535" = 3,
               "SRR14800536" = 3,
               "SRR14800540" = 4,
               "SRR14800541" = 3,
               "SRR14800543" = 3,
               "SRR17960481" = 3,
               "SRR17960484" = 3
             )),

  tar_target(volcano_large_in_segment,
             make_volcano_diffex_clones(large_in_segment_diffex_clones_for_each_cluster,
                                    "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
                                    large_in_segment_diffex_clones,
                                    "results/diffex_bw_clones_large_in_segment.pdf"),
  ),

  tar_target(table_large_in_segment_diffex_clones,
             tabulate_diffex_clones(large_in_segment_diffex_clones_for_each_cluster,
                                    "results/diffex_bw_clones_per_cluster_large_in_segment.xlsx",
                                    "results/diffex_bw_clones_per_cluster_large_in_segment_by_chr.xlsx",
                                    large_in_segment_diffex_clones,
                                    "results/diffex_bw_clones_large_in_segment.xlsx",
                                    "results/diffex_bw_clones_large_in_segment_by_chr.xlsx"),
  ),

  # out of segment ------------------------------
  tar_target(large_out_of_segment_diffex_clones,
             find_diffex_clones(large_done_files, large_clone_comparisons, location = "out_of_segment"),
             pattern = map(large_done_files),
             iteration = "list"
  ),

  tar_target(
    large_out_of_segment_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(large_done_files, large_clone_comparisons, cluster_dictionary, location = "out_of_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),

  tar_target(volcano_large_out_of_segment,
             make_volcano_diffex_clones(large_out_of_segment_diffex_clones_for_each_cluster,
                                        "results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf",
                                        large_out_of_segment_diffex_clones,
                                        "results/diffex_bw_clones_large_out_of_segment.pdf"),
  ),


  tar_target(table_large_out_of_segment_diffex_clones,
             tabulate_diffex_clones(large_out_of_segment_diffex_clones_for_each_cluster,
                                    "results/diffex_bw_clones_per_cluster_large_out_of_segment.xlsx",
                                    "results/diffex_bw_clones_per_cluster_large_out_of_segment_by_chr.xlsx",
                                    large_out_of_segment_diffex_clones,
                                    "results/diffex_bw_clones_large_out_of_segment.xlsx",
                                    "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx"),
  ),


  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(large_in_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_in_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  #
  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(large_out_of_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_out_of_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(large_in_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_in_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  #

  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(large_out_of_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_out_of_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  tar_target(
    large_diffex_bw_clusters_for_each_clone,
    find_diffex_bw_clusters_for_each_clone(large_done_files, cluster_dictionary),
    pattern = map(large_done_files),
    iteration = "list"
  ),

#
# tar_target(
#   large_collated_heatmaps,
#   qpdf::pdf_combine(large_heatmaps, "results/numbat_sridhar_large/large_heatmaps.pdf"),
# ),
#
tar_target(large_numbat_expression, retrieve_numbat_plot_type(large_numbat_pdfs, "exp_roll_clust.pdf")),

# tar_target(large_clone_dist, retrieve_numbat_plot_type(large_plot_files, "exp_roll_clust.pdf")),

# tar_target(large_numbat_expression_pdf, qpdf::pdf_combine(large_numbat_expression, "results/numbat_sridhar_large/large_numbat_expression.pdf")),
#
# tar_target(large_numbat_heatmap_files, retrieve_numbat_plot_type(large_numbat_pdfs, "panel_2.pdf")),
#
# tar_target(large_numbat_heatmap_pdf, qpdf::pdf_combine(large_numbat_heatmap_files, "results/numbat_sridhar_large/large_numbat_heatmaps.pdf")),
#
# tar_target(large_numbat_bulk_clones_final, retrieve_numbat_plot_type(large_numbat_pdfs, "bulk_clones_final.pdf")),
#
# tar_target(large_numbat_bulk_clones_final_pdf, qpdf::pdf_combine(large_numbat_bulk_clones_final, "results/numbat_sridhar_large/large_bulk_clones_final.pdf")),
#
# # large_numbat_sample_pdfs ------------------------------
tar_target(
  large_numbat_sample_pdfs,
  reroute_done_to_results_pdf(large_done_files, "_large"),
  pattern = map(large_done_files),
  iteration = "list",
  format = "file"
),
#
# tar_target(
#   large_cluster_markers,
#   collect_all_markers(large_done_files, "results/clusters.xlsx")
# )
#
# # tar_target(
# #   large_marker_gene_featureplots,
# #   plot_putative_marker_across_samples(interesting_genes, large_done_files, plot_type = FeaturePlot),
# # ),
#
# # tar_target(
# #   large_marker_gene_vlnplots_by_cluster,
# #   plot_putative_marker_across_samples(interesting_genes, large_done_files, VlnPlot, group_by = "abbreviation", cluster_dictionary),
# # ),
# #
# # tar_target(
# #   large_marker_gene_vlnplots_by_clone,
# #   plot_putative_marker_across_samples(interesting_genes, large_done_files, VlnPlot, group_by = "clone_opt", cluster_dictionary),
# # ),
# end of plan------------------------------
)

# notes 2023-03-07 ------------------------------

  # yes ------------------------------
  # wu SRR13884242 (RB2 rep 1): likely 16q GT; f31 sample; cluster 1 (prolif.) has higher proportion of GT 2; can attribute this proliferation to acquisition of 16q- in GT 2
  # wu SRR13884243 (RB2 rep 2): likely 16q GT; biological replication of SRR13884242; cluster 4 is an analogue of SRR13884242 c3

  # wu SRR13884247: likely 16q GT; no GSEA diffex output; cluster 5 is interesting; c5 marker DLK1 in amicroRNA cluster with MEG3
  # wu SRR13884249: possible 1q/2p/16q GT; evidence that GT 1 (only 16q-) does not contribute to proliferating clusters c2 and c4 (TOP2A high); split by phase
  # yang SRR14800534: GT1 lacks 16q-; does not contribute to proliferating clusters c2 and c4; GT 2 lacks 1q--can't identify tx distinction b/w gt2/3
  # yang SRR14800535: 16q GT; GT 1 decreased in c1; c1 high TOP2A
  # yang SRR14800536: possible 16q GT; GT1 decreased in c2/3 high G2M markers
  # yang SRR14800540: clear 16q GT; very complicated tumor; three GTs with SCNAs; each is proliferating to a greater degree; confused about GT5 with high c4; markers are C1QA/B, CD74, HLA genes
  # yang SRR14800541: clear 1q 6p and 16q GTs; GT1 not proliferating--no contribution to c2 with G2M markers
  # yang SRR14800543: possible minor 16q/1q GT; can identify dual or individual contribution of 1q and 16q with 13q CNLOH; strange MYCN marker of clusters
  # field SRR17960481: 6p interesting; cluster 2 and 4; cluster 2 has mito genes; not promising; likely clonal 6p with PRs and stressed cells composing GT 1
  # field SRR17960484: cluster 1 enriched for wt GT 1; also has high Xist expression

  # maybe ------------------------------
  # SRR13884240: possible 2p GT
  # SRR13884241: possible 1q GT
  # SRR13884244: possible 1q GT
  # SRR13884245: possible 1q GT
  # SRR13884246: possible 16q GT
  # SRR17960480: possible minor 16q- GT
  # SRR14800539: possible 16q GT

  # no ------------------------------
  # SRR14800537: possible 16q GT
  # SRR17960482: too complicated
  # SRR13884248: clear 6p (missing possible 2p in expression)
  # SRR17960483: cluster 6 maybe interesting
