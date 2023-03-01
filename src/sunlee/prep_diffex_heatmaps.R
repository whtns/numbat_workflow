#!/usr/bin/Rscript

# now let's spice up the dendrograms a bit:
library(dendextend)
library(tidyverse)
library(data.table)


census_matrix = read.table("~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/sunhye_census_matrix_20170407.csv", sep = "\t", header = TRUE)

annotation = read.table("~/single_cell_pipeline/scde_input/shl_0407_w_centroids_cell_info.csv", sep = "\t", header = TRUE)


diffex_csvs <- list.files("/home/skevin/single_cell_tools/p_val_filtered_scde/", full.name = TRUE)

cluster_names <- gsub(".csv", "", basename(diffex_csvs))

cluster_diffex <- lapply(diffex_csvs, read.table, sep = ",", header = TRUE, row.names =1)
names(cluster_diffex) <- cluster_names

cluster_diffex_upr <- purrr::map(cluster_diffex, dplyr::top_n, 20, mle)
cluster_diffex_dnr <- purrr::map(cluster_diffex, dplyr::top_n, -20, mle)

cluster_diffex <- mapply(rbind, cluster_diffex_upr, cluster_diffex_dnr, SIMPLIFY=FALSE)


cluster_733_234 <- dplyr::bind_rows(cluster_diffex[c(seq(1,6))])
cluster_733_345 <- dplyr::bind_rows(cluster_diffex[c(seq(7,12))])
cluster_737 <- dplyr::bind_rows(cluster_diffex[c(seq(13,18))])
cluster_ctrl <- dplyr::bind_rows(cluster_diffex[c(seq(19,24))])

cluster_diffex <- list("733_234" = cluster_733_234, "733_345" = cluster_733_345, "737" = cluster_737, "ctrl" = cluster_ctrl)

new_cluster_names <- c("733_234", "733_345", "737", "ctrl")

filter_cms <- lapply(cluster_diffex, function(x) {census_matrix[x[,"ensembl_transcript_id"],]})

treatment_levels = levels(annotation$treatment_group)

annotations <- lapply(treatment_levels, function(x) annotation[annotation["treatment_group"]==x,])

annotations <- c(annotations[1], annotations[1], annotations[2], annotations[3])

cluster_cols <- names(annotation)[c(9,10,7,8)]

cluster_sets0 <- purrr::map2(annotations, cluster_cols, function(x,y) x[!is.na(x[[y]]),])
cluster_dictdfs <- purrr::map2(cluster_sets0, cluster_cols, function(x,y) dplyr::select(x, c("Sample_ID", "cluster" = y)))


# cluster_cols <- quos(cluster_cols)


# score[order(sex, y, x),]
# 
# cluster_sets <- purrr::map2(annotations, cluster_cols, function(x,y) arrange(x, !!! y))

annotations[[1]] <- arrange(cluster_sets[[1]], centroid_733_234)
annotations[[2]] <- arrange(cluster_sets[[2]], centroid_733_345)
annotations[[3]] <- arrange(cluster_sets[[3]], centroid_737)
annotations[[4]] <- arrange(cluster_sets[[4]], centroid_ctrl)


#combine diffex from same analysis different cluster




# filter_cms <- purrr::map2(filter_cms, annotations, function(x,y) x[,colnames(x) %in% y[["Sample_ID"]]])

filter_cms <- purrr::map2(filter_cms, annotations, function(x,y) dplyr::select(x, as.character(y[["Sample_ID"]])))

filter_cms_trs <- filter_cms
filter_cms_genes <- filter_cms

purrr::map2(filter_cms_genes, cluster_diffex, function(x,y) {row.names(x) <- make.names(y[["hgnc_symbol"]], unique = TRUE)
                                        return(x)})

genes_names <- paste0(new_cluster_names, "_genes_filtered.csv")
trs_names <- paste0(new_cluster_names, "_transcripts_filtered.csv")

# filter_cms <- lapply(filter_cms, rownames_to_column, var = "ensembl_transcript_id")
purrr::map2(filter_cms_genes, genes_names, write.table, sep=",")
purrr::map2(filter_cms_trs, trs_names, write.table, sep=",")

filter_cms_genes <- lapply(genes_names, read.csv)

# write.table(filt_census_matrix, "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/scde_input/diffex_by_trs_clusters_1_4/heatmap_input.csv")

# filt_cms  <- lapply(filter_cms, function(x) data.matrix(x))


test <- filt_cms[[1]]

col_dend  <- test %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k = 2) %>% set("branches_lwd", c(1,2)) %>%
  ladderize
#    rotate_DendSer(ser_weight = dist(t(x)))

heatmaply(percentize(x), Colv = col_dend)