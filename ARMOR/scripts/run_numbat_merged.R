#!/usr/bin/env Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

library(numbat)
library(Seurat)
library(readr)
library(magrittr)
conflicted::conflict_prefer("rowSums", "Matrix")

allele_df = read_tsv("../output/numbat/merged/merged_allele_counts.tsv.gz") %>% 
  tidyr::drop_na(CHROM) %>% # drop X chromosome values %>% 
	dplyr::mutate(CHROM = as.factor(CHROM)) %>% 
  identity()

count_mat <- readRDS("../output/numbat/count_mat.rds")

out_dir = "../output/numbat/merged/"
ncores = "6"

# mygtf_hg38 <- readRDS("~/tmp/gtf_transcript.rds")

# debug(numbat::run_numbat)

# run
out = run_numbat(
	count_mat, # gene x cell integer UMI count matrix done
	ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix done
	allele_df, # allele dataframe generated by pileup_and_phase script done
	gtf_hg38, # provided upon loading the package done
	genetic_map_hg38, # provided upon loading the package done
	min_cells = 10,
	t = 1e-3,
	max_iter = 2,
	max_entropy = 0.60,
	min_LLR = 50,
	init_k = 3,
	ncores = as.integer(ncores),
	plot = TRUE,
	multi_allelic = FALSE,
	common_diploid = TRUE,
	out_dir = out_dir,
	tau = 0.3)
