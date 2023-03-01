#!/usr/bin/Rscript

#load required libraries and functions
#=====================================
library(optparse)

# load( "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_0407_deseq_input.rda")
# list2env(shl, globalenv())

#SHL 20171031
# default_infile <- "~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/FACS_20171031_sunlee_H_sapiens_summarized_experiment.rds"
# default_cell_info = "~/single_cell_tools/FACS_0407_1031_SHL_input_files/cells_sets_1_2.csv"
# default_out = "~/monocle_test.pdf"

#SHL 20170407
default_infile <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/FACS_20170407_sunlee_H_sapiens_summarized_experiment.rds"
default_ptime <- "~/tmp/5sh733_shCtrl_PT"
default_tdf_location <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/igv_tdf_files/"

#DS 20171031_20170407
# default_infile <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_summarized_experiment.rds"
# default_cell_info = "~/single_cell_tools/dshayler_input/2_seq_3dformat_050418.csv"
# default_plot_settings = "~/single_cell_tools/dshayler_input/030618_3d_PCA_No_Bad_Reads_Color_by_age.txt"
# default_out = "~/ds_test.pdf"

# default input files -----------------------------------------------------

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=default_infile,
              help="summarized experiment input filename [default= %default]", metavar="character"),
  make_option(c("-p", "--pseudotime_file"), type="character", default=default_ptime,
              help="pseudotime file for ordering cells [default= %default]", metavar="character"),
  make_option(c("-t", "--tdf_location"), type="character", default=NA,
              help="directory that holds all tdf files [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  missing_args <- paste(names(which(sapply(opt, is.na))), collapse = " ")
  stop(paste0("Please provide all necessary arguments. Missing: ", missing_args), call.=FALSE)
}

# load required libraries -------------------------------------------------
print("loading required libraries")
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))

# create attributes file --------------------------------------------------
sce <- readRDS(opt$infile)

cd <- SummarizedExperiment::colData(sce)

ptime <- read.table(opt$pseudotime_file)
colnames(ptime) <- c("sample_id", "pseudotime")

tdf_files <- list.files(opt$tdf_location, pattern = "*.tdf")
tdf_prefixes <- paste0("X", gsub("_.*", "", tdf_files))

tdf_df <- data.frame("Array" = tdf_files, "sample_id" = tdf_prefixes)

igv_attributes <- as.data.frame(cd) %>% 
  right_join(ptime, by = "sample_id")  %>% 
  dplyr::arrange(treatment_group, pseudotime) %>% 
  left_join(tdf_df, by = "sample_id") %>% 
  dplyr::mutate(Order = 1:length(sample_id)) %>% 
  dplyr::select(Array, Order, treatment_group, pseudotime)

out_txt <- paste0(opt$tdf_location, "/igv_attributes.txt")

print(paste0("saving attributes file as: ", out_txt))
write.table(igv_attributes, out_txt, row.names = FALSE, sep = "\t", quote = F)



