#!/usr/bin/Rscript


# load optparse -----------------------------------------------------------
library(optparse)

# handle command line options ---------------------------------------------

#default input files -----------------------------------------------------
# default_expression_files = system.file("extdata", "expression_matrix_sample.csv", package="cataract")
# default_cell_info = system.file("extdata", "coldata_sample.csv", package="cataract")
# default_out = "./scde"

default_expression_files = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/Sunhye_stringtie.tpm.csv:~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/transcripts.tpm.csv"
default_cell_info = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/Sunhye_cell_division_day_treatment.csv:~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20171031_sunlee_H_sapiens/FACS_20171031_sunlee_sample_sheet.csv"
default_out = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_sunlee_FACS_20171031_sunlee_singlecellexperiment2.rds"
default_tags = "exp1_,exp2_"

# default_expression_files = "~/single_cell_pipeline/output/dshayler/20170407_20171031_combined_stringtie.tpm.csv,~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/transcripts.tpm.csv" 
# default_cell_info = "~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/2_seq_dshayler/2_seq_meta_111517.csv,~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20171031_sunlee_H_sapiens/FACS_20171031_sunlee_sample_sheet-tab.csv"
# default_out = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_dshayler_FACS_20171031_sunlee_singlecellexperiment2.rds"
# default_tags = "DS_,SHL_"

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-e", "--gene_expression_files"), type="character", default=default_expression_files,
              help="colon separated list of gene expression input filenames [default= %default]", metavar="character"),
  make_option(c("-c", "--cell_info"), type="character", default=default_cell_info,
              help="colon separated list of cell metadata filenames [default= %default]", metavar="character"),
  make_option(c("-t", "--exp_tags"), type="character", default=NA,
              help="colon separated list of experiment tags (ex. DS_:SHL_) [default= %default]", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=default_out,
              help="output_file for merged singlecellexperiment [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

# load required libraries -------------------------------------------------
library(reshape2)
require(gtools)
library(tidyverse)
library(SingleCellExperiment)


# load required data ------------------------------------------------------

EXPRESSION_PATHS <- unlist(strsplit(opt$gene_expression_files,"[,:]"))
GROUP_PATHS = unlist(strsplit(opt$cell_info,"[,:]"))
exp_tags = unlist(strsplit(opt$exp_tags,"[,:]"))
outfile = opt$outfile

safe_read <- function(sep_file){
  L <- readLines(sep_file, n = 1)
  if (grepl("\t", L)) read.table(sep_file, sep = "\t", header = TRUE) else read.csv(sep_file)
}

counts <- purrr::map(EXPRESSION_PATHS, safe_read)
colData <- purrr::map(GROUP_PATHS, safe_read)

#check if dimensions of counts and metadata match!
test_d_m <- purrr::map2(counts, colData, function(x,y) all.equal(ncol(x), nrow(y), tolerance = 1))
if(!isTRUE(all(test_d_m))){
  print("data and metadata have unequal number of samples. Check input files")
  quit()
  
}

# construct a merged SingleCellExperiment ---------------------------------
make_sce <- function(counts, colData, exp_tags){

  check_rownames <- function(counts){
    if(grepl("ENST.*", rownames(counts))){
      rownames(counts) <- rownames(counts)
    } else{
      rownames(counts) <- counts[,1]
      counts[,1] <- NULL
    } 
    return(counts)
  }
  
  counts <- lapply(counts, check_rownames)
  
  data_append_exp_tag <- function(df, exp_tag){
    names(df) <- paste0(exp_tag, names(df))
    return(df)
  }
  
  meta_append_exp_tag <- function(df, exp_tag){
    names(df)[1] <- "sample_id"
    df <- dplyr::mutate(df, sample_id = paste0(exp_tag, sample_id))
    return(df)
  }
  
  names(counts) <- exp_tags
  names(colData) <- exp_tags
  counts <- purrr::map2(counts, exp_tags, data_append_exp_tag)
  colData <- purrr::map2(colData, exp_tags, meta_append_exp_tag)
  
  if (dim(counts[[1]]) > dim(counts[[2]])){
    counts <- cbind(counts[[2]], counts[[1]][match(rownames(counts[[2]]), rownames(counts[[1]])),])
  } else if (dim(counts[[1]] < dim(counts[[2]]))){
    counts <- cbind(counts[[1]], counts[[2]][match(rownames(counts[[1]]), rownames(counts[[2]])),])
  }
  colData <- dplyr::bind_rows(colData)
  colData <- colData[order(colData$sample_id),]
  # rownames(colData) <- colData["sample_id"]
  counts <- counts[,order(colnames(counts))]
  counts <- as.matrix(counts)
  census_experiment <- SingleCellExperiment(assays=list(counts=counts), colData=colData, rowData=rownames(counts))
  return(census_experiment)
}

sce <- make_sce(counts, colData, exp_tags)

saveRDS(sce, outfile)
