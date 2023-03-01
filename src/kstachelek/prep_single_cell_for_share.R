#!/usr/bin/Rscript

library(dplyr)
library(SingleCellExperiment)




raw_counts_f <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/transcript_count_matrix.csv"

fpkm_f <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/Sunhye_stringtie.fpkm.csv"

census_f <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/sunhye_census_matrix_20170407.csv"

raw_counts <- read.csv(raw_counts_f)
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]

fpkm <- read.table(fpkm_f, sep="\t", header = TRUE)
rownames(fpkm) <- fpkm[,1]
fpkm <- fpkm[,-1]

census <- read.table(census_f, sep="\t", header = TRUE)

test <- lapply(c(raw_counts_f, fpkm_f, census_f), cataract::safe_read)

meta_f <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/Sunhye_cell_division_day_treatment.csv"



meta <- read.csv(meta_f)

table_names <- c("raw_counts", "fpkm", "normalized", "metadata")

outfiles <- setNames(list(raw_counts, fpkm, census, meta), table_names)

file_names <- paste0("~/sunlee_cobrinik_", table_names, ".csv")

purrr::map2(outfiles, file_names, write.csv)


