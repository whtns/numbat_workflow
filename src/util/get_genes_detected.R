#!/usr/bin/Rscript

library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)

#SHL 20170407
# transcript_counts <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/transcript_count_matrix.csv"
# gene_counts <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/gene_count_matrix.csv"

#SHL 20171031
# transcript_counts <- "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/transcript_count_matrix.csv"
# gene_counts <- "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/gene_count_matrix.csv"

meetsCondition = function(mat, gtet = 0, nCells=0){
  condition_mat = ifelse(mat>=gtet,1,0)
  meets_condition = apply(condition_mat,1,sum) >= nCells
  return(meets_condition)
}


#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-c", "--countsfile"), type="character", default=NA,
              help="gene expression input filename [default= %default]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NA,
              help="output directory name [default= %default]", metavar="character"),
  make_option(c("-t", "--threshold"), type="integer", default=NA, 
              help="the minimum number of cells in which a gene must be expressed to be called detected", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

counts <- data.frame(fread(opt$countsfile, header = TRUE), row.names = 1)

counts_over_threshold <- counts[meetsCondition(counts, 1, opt$threshold),]


dir.create(opt$outdir, showWarnings = FALSE)

genes_detected <- data.frame("genes_detected" = colSums(counts != 0)) %>% 
  tibble::rownames_to_column("cell_id")
colnames(genes_detected) <- c("cell_id", "genes_detected")
write.csv(genes_detected, paste0(opt$outdir, "genes_detected_in_at_least_one_cell.csv"), row.names = FALSE)

if(!opt$threshold == 1){
  genes_detected_min_cells <- data.frame("genes_detected" = colSums(counts_over_threshold != 0)) %>% 
    tibble::rownames_to_column("cell_id")
  colnames(genes_detected) <- c("cell_id", paste0("genes_detected", "_in_", as.character(opt$threshold), "_cells"))
  write.csv(genes_detected_min_cells, paste0(opt$outdir, "genes_detected_in_at_least_", as.character(opt$threshold), "_cells.csv"), row.names = FALSE)
}
