#!/usr/local/bin/Rscript

default_cell_info <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/cell_info.csv"

outdir <- "/home/sunlee/igv_tdf_files/"

cell_info <- read.csv(default_cell_info)

cell_info <- cell_info[!is.na(cell_info$scheme1),]

cluster_dfs <- split(cell_info, cell_info$scheme1)

cluster_ids <- lapply(cluster_dfs, function(x) x[["Sample_ID"]] = paste0(gsub("X", "/", x[["Sample_ID"]]), "_"))

patterns <- lapply(cluster_ids, paste, sep="", collapse="|")

remove_dup_bam_files <- list.files("~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/", pattern = ".*duplicates.bam$", recursive  = TRUE, full.names = TRUE)

cluster_files <- lapply(patterns, function(x) remove_dup_bam_files[grepl(x, remove_dup_bam_files)])

igv_cluster_files <- paste0(outdir, names(cluster_files), ".bam.list")


purrr::map2(cluster_files, igv_cluster_files, write.table, row.names = FALSE, col.names = FALSE, quote = FALSE)
