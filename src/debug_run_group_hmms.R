#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)

bulk_subtrees_2_4242 <- read_tsv("output/numbat_sridhar/SRR13884242/bulk_subtrees_2.tsv.gz") %>%
  dplyr::mutate(CHROM = factor(CHROM))

bulk_clones_final_4242 <- read_tsv("output/numbat_sridhar/SRR13884242/bulk_clones_final.tsv.gz") %>%
  dplyr::mutate(CHROM = factor(CHROM))

# debug(numbat:::run_group_hmms)

bulk_subtrees <- numbat:::run_group_hmms(bulk_subtrees_2_4242)

bulk_clones <- numbat:::run_group_hmms(bulk_clones_final_4242)

# test0 <-
#   bulks %>%
#   dplyr::filter(CHROM == "14")

p_subtress = numbat:::plot_bulks(bulk_subtrees, min_LLR = 2, use_pos = TRUE, genome = 'hg38')

p_clones_final = numbat:::plot_bulks(bulk_clones, min_LLR = 2, use_pos = TRUE, genome = 'hg38')
