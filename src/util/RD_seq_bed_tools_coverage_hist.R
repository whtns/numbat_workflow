#!/bin/Rscript

library(ggplot2)
library(dplyr)

unfiltered_path <- "~/single_cell_pipeline/output/GT_20171101_20171016_merged_SHL_H_sapiens_output/test.txt"
intronic_filtered_path <- "~/single_cell_pipeline/output/GT_20171101_20171016_merged_SHL_H_sapiens_output/test_unique_intronic_2.txt" 
depth_under_3_path <- "~/single_cell_pipeline/output/GT_20171101_20171016_merged_SHL_H_sapiens_output/test_depth_under_3_2.txt"

prep_df <- function(hist_path, hist_plt_path){
  browser()
  hist <- read.table(hist_path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  header <- c("chrom", "start", "end", "read_id", "score", "strand", "thick_start", "thick_end", "rgb", "blckcount", "blcksizes", "blckstarts", "depth", "bases_at_depth", "size", "fraction_covered")
  names(hist) <- header
  
  hist <- mutate(hist, exon_overlap = (bases_at_depth/100))
  
  hist_plt <- ggplot(hist, aes(x=fraction_covered, y=depth)) + geom_point()
  ggsave(hist_plt, hist_plt_path)
  return(hist)
}

unfiltered <- prep_df(unfiltered_path, "~/single_cell_pipeline/results/kstachelek/lab_meeting_20171215/unfiltered_depth_by_coverage.pdf")

intronic <- prep_df(intronic_filtered_path, "~/single_cell_pipeline/results/kstachelek/lab_meeting_20171215/intronic_filtered_depth_by_coverage.pdf" )

depth_under_3 <- prep_df(depth_under_3_path, "~/single_cell_pipeline/results/kstachelek/lab_meeting_20171215/depth_under_3_depth_by_coverage.pdf" )



unfiltered_plt <- ggplot(unfiltered, aes(x=fraction_covered, y=depth)) + geom_point()  

dir.create("~/single_cell_pipeline/results/kstachelek/lab_meeting_20171215", recursive = TRUE)
ggsave(unfiltered, "~/single_cell_pipeline/results/kstachelek/lab_meeting_20171215/unfiltered_depth_by_coverage.pdf")
