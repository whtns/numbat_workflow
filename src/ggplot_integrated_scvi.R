#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(patchwork)
library(ggalluvial)


order_df_by_diversity <- function(df, group1, group2){

  test0 <-
    df %>%
    group_by(.data[[group1]], .data[[group2]]) %>%
    dplyr::summarize(values = n()) %>%
    group_by(.data[[group2]]) %>%
    summarize(
      shannon = vegan::diversity(values, index = "shannon")
    ) %>%
    arrange(shannon) %>%
    identity()

  res_df <-
    df %>%
    left_join(test0, by = c(group2)) %>%
    arrange(shannon) %>%
    dplyr::mutate({{group2}} := factor(.data[[group2]], levels = unique(.data[[group2]]))) %>%
    identity()

  return(res_df)

}

plot_metadata_category_patchwork_barplot <- function(metadata_path = "results/adata_all_metadata.csv", plot_path = "results/barplot_patchworks.pdf"){

  combined_adata_obs <-
    metadata_path %>%
    read_csv() %>%
    dplyr::rename(cell = `...1`) %>%
    dplyr::mutate(cluster_short = str_replace(cluster, ".*_", "")) %>%
    dplyr::mutate(cluster_short = dplyr::case_when(cluster_short == "ARL1IP1" ~ "ARL6IP1",
                                                   TRUE ~ cluster_short)) %>%
    dplyr::filter(!is.na(cluster_short))


  leiden_v_sample_plot <-
    combined_adata_obs %>%
    order_df_by_diversity("sample_id", "leiden") %>%
    dplyr::mutate(leiden = factor(leiden)) %>%
    ggplot(aes(x = leiden, fill = sample_id)) +
    geom_bar(position = "fill") +
    coord_flip()

  cluster_v_leiden_plot <-
    combined_adata_obs %>%
    order_df_by_diversity("leiden", "cluster_short") %>%
    dplyr::mutate(leiden = factor(leiden)) %>%
    ggplot(aes(x = cluster_short, fill = leiden)) +
    geom_bar(position = "fill") +
    coord_flip()

  leiden_v_cluster_plot <-
    combined_adata_obs %>%
    order_df_by_diversity("cluster_short", "leiden") %>%
    ggplot(aes(x = leiden, fill = cluster_short)) +
    geom_bar(position = "fill") +
    coord_flip()

  cluster_v_sample_plot <-
    combined_adata_obs %>%
    order_df_by_diversity("sample_id" ,"cluster_short") %>%
    ggplot(aes(x = cluster_short, fill = sample_id)) +
    geom_bar(position = "fill") +
    coord_flip()

  cluster_v_sample_plot + leiden_v_sample_plot + cluster_v_leiden_plot + leiden_v_cluster_plot

  ggsave(plot_path, width = 12, height = 8)
}

plot_metadata_category_patchwork_barplot(metadata_path = "results/adata_all_metadata.csv", plot_path = "results/cluster_v_sample_v_leiden_barplot_all.pdf")

plot_metadata_category_patchwork_barplot(metadata_path = "results/adata_filtered_metadata.csv", plot_path = "results/cluster_v_sample_v_leiden_barplot_filtered.pdf")

# alluvial plots ------------------------------

plot_metadata_category_patchwork_alluvial <- function(metadata_path = "results/adata_all_metadata.csv", plot_path = "results/alluvial_patchworks.pdf", cluster_order = c("MYCN", "HSP", "MEG3", "MALAT1", "HIST1H4C", "PCLAF", "MARCKSL1", "TFF1", "ARL6IP1", "cone", "RACK1", "G2M")){

  # browser()

  cluster_v_leiden_plot <-
    metadata_path %>%
    read_csv() %>%
    dplyr::rename(cell = `...1`) %>%
    dplyr::mutate(cluster_short = str_replace(cluster, ".*_", "")) %>%
    dplyr::mutate(cluster_short = dplyr::case_when(cluster_short == "ARL1IP1" ~ "ARL6IP1",
                                                   TRUE ~ cluster_short)) %>%
    dplyr::filter(!is.na(cluster_short)) %>%
    # dplyr::mutate(cluster_short = factor(cluster_short, levels = cluster_order)) %>%
    dplyr::mutate(leiden = factor(leiden)) %>%
    group_by(cluster_short, leiden) %>%
    dplyr::summarize(Freq = n()) %>%
    ggplot(aes(axis1 = cluster_short, axis2 = leiden, y = Freq)) +
    scale_x_discrete(limits = c("cluster_short", "leiden"), expand = c(.2, .05)) +
    xlab("cluster") +
    geom_alluvium(aes(fill = cluster_short),
                  aes.bind = TRUE,
                  decreasing = TRUE) +
    geom_stratum(decreasing = TRUE) +
    geom_text(stat = "stratum", decreasing = TRUE, aes(label = after_stat(stratum))) +
    theme_minimal()

  leiden_v_cluster_plot <-
    metadata_path %>%
    read_csv() %>%
    dplyr::rename(cell = `...1`) %>%
    dplyr::mutate(cluster_short = str_replace(cluster, ".*_", "")) %>%
    dplyr::mutate(cluster_short = dplyr::case_when(cluster_short == "ARL1IP1" ~ "ARL6IP1",
                                                   TRUE ~ cluster_short)) %>%
    dplyr::filter(!is.na(cluster_short)) %>%
    # dplyr::mutate(cluster_short = factor(cluster_short, levels = cluster_order)) %>%
    dplyr::mutate(leiden = factor(leiden)) %>%
    group_by(cluster_short, leiden) %>%
    dplyr::summarize(Freq = n()) %>%
    ggplot(aes(axis1 = leiden, axis2 = cluster_short, y = Freq)) +
    scale_x_discrete(limits = c("cluster_short", "leiden"), expand = c(.2, .05)) +
    xlab("cluster") +
    geom_alluvium(aes(fill = leiden),
                  aes.bind = TRUE) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()

  # leiden_v_sample_plot <-
  #   metadata_path %>%
  #   read_csv() %>%
  #   dplyr::rename(cell = `...1`) %>%
  #   dplyr::mutate(cluster_short = str_replace(cluster, ".*_", "")) %>%
  #   dplyr::mutate(cluster_short = dplyr::case_when(cluster_short == "ARL1IP1" ~ "ARL6IP1",
  #                                                  TRUE ~ cluster_short)) %>%
  #   dplyr::filter(!is.na(cluster_short)) %>%
  #   group_by(sample_id, leiden) %>%
  #   dplyr::summarize(Freq = n()) %>%
  #   ggplot(aes(axis1 = sample_id, axis2 = leiden, y = Freq)) +
  #   scale_x_discrete(limits = c("sample_id", "leiden"), expand = c(.2, .05)) +
  #   xlab("cluster") +
  #   geom_alluvium(aes(fill = sample_id)) +
  #   geom_stratum() +
  #   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  #   theme_minimal() +
  #   ggtitle("clusters")
  #
  # cluster_v_sample_plot <-
  #   metadata_path %>%
  #   read_csv() %>%
  #   dplyr::rename(cell = `...1`) %>%
  #   dplyr::mutate(cluster_short = str_replace(cluster, ".*_", "")) %>%
  #   dplyr::mutate(cluster_short = dplyr::case_when(cluster_short == "ARL1IP1" ~ "ARL6IP1",
  #                                                  TRUE ~ cluster_short)) %>%
  #   dplyr::filter(!is.na(cluster_short)) %>%
  #   group_by(sample_id, cluster_short) %>%
  #   dplyr::summarize(Freq = n()) %>%
  #   ggplot(aes(axis1 = cluster_short, axis2 = sample_id, y = Freq)) +
  #   scale_x_discrete(limits = c("cluster_short", "sample_id"), expand = c(.2, .05)) +
  #   xlab("cluster") +
  #   geom_alluvium(aes(fill = sample_id)) +
  #   geom_stratum() +
  #   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  #   theme_minimal() +
  #   ggtitle("clusters")

  # cluster_v_sample_plot +
  #   leiden_v_sample_plot +
    # cluster_v_leiden_plot

  cluster_v_leiden_plot + leiden_v_cluster_plot

  ggsave(plot_path, width = 16, height = 12)
}

plot_metadata_category_patchwork_alluvial(metadata_path = "results/adata_all_metadata.csv", plot_path = "results/cluster_v_leiden_alluvial_all.pdf")

plot_metadata_category_patchwork_alluvial(metadata_path = "results/adata_filtered_metadata.csv", plot_path = "results/cluster_v_leiden_alluvial_filtered.pdf")


  # grouped subsets ------------------------------

# plot_metadata_category_patchwork_barplot(metadata_path = "results/adata_16q_plus_1q_metadata.csv", plot_path = "results/cluster_v_sample_v_leiden_barplot_16q_plus_1q.pdf")
#
# plot_metadata_category_patchwork_barplot(metadata_path = "results/adata_1q_6p_metadata.csv", plot_path = "results/cluster_v_sample_v_leiden_barplot_1q_6p.pdf")
#
# plot_metadata_category_patchwork_barplot(metadata_path = "results/adata_16q_metadata.csv", plot_path = "results/cluster_v_sample_v_leiden_barplot_16q.pdf")
#
# plot_metadata_category_patchwork_barplot(metadata_path = "results/adata_1q_16q_metadata.csv", plot_path = "results/cluster_v_sample_v_leiden_barplot_1q_16q.pdf")

# plot_metadata_category_patchwork_alluvial(metadata_path = "results/adata_16q_plus_1q_metadata.csv", plot_path = "results/cluster_v_leiden_alluvial_16q_plus_1q.pdf")
#
# plot_metadata_category_patchwork_alluvial(metadata_path = "results/adata_1q_6p_metadata.csv", plot_path = "results/cluster_v_leiden_alluvial_1q_6p.pdf")
#
# plot_metadata_category_patchwork_alluvial(metadata_path = "results/adata_16q_metadata.csv", plot_path = "results/cluster_v_leiden_alluvial_16q.pdf")
#
# plot_metadata_category_patchwork_alluvial(metadata_path = "results/adata_1q_16q_metadata.csv", plot_path = "results/cluster_v_leiden_alluvial_1q_16q.pdf")

# tables------------------------------

# browseURL("results/cluster_v_sample_v_leiden.pdf")
#
#
# cluster_table <-
#   combined_adata_obs %>%
#   group_by(leiden, cluster_short) %>%
#   summarize(n_cluster_cells = n()) %>%
#   group_by(leiden) %>%
#   mutate(per =  100 *n_cluster_cells/sum(n_cluster_cells)) %>%
#   # ungroup %>%
#   dplyr::arrange(leiden, desc(per)) %>%
#     identity()
#
#
# sample_table <-
#   combined_adata_obs %>%
#   group_by(leiden, sample_id) %>%
#   summarize(n_cluster_cells = n()) %>%
#   group_by(leiden) %>%
#   mutate(per =  100 *n_cluster_cells/sum(n_cluster_cells)) %>%
#   # ungroup %>%
#   dplyr::arrange(leiden, desc(per)) %>%
#   identity()
