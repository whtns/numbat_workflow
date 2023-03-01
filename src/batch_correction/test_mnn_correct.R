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
  make_option(c("-t", "--exp_tags"), type="character", default=default_tags,
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
library(scran)

mt_settings <- convert_mt_setting(opt$cellset, opt$plot_settings)
# load required functions -------------------------------------------------

## parse martin metadata
find_remove_cells <- function(plot_settings, annotation){
  
  test <- readLines(plot_settings)
  
  if (!grepl('remove', test)){
    return(NULL)
  }
  
  vecs <- list()
  mtnames <- c()
  for (i in test){
    if (!grepl("#", i) & grepl('remove', i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[2]]
      vecs <- append(vecs, lbline)
      mtnames <- append(mtnames, d)    
    }
  }
  
  pfx <- tolower(gsub("_.*", "", mtnames))
  valid_pfx  <- which(pfx %in% colnames(annotation))
  
  if (length(valid_pfx) == 0){
    return(NULL)
  }
  
  pfx <- pfx[valid_pfx]
  
  
  
  sfx <- gsub(".*_", "", mtnames[valid_pfx])
  
  remove_cells <- purrr::map2(pfx, sfx , function(x, y) annotation[annotation[tolower(x)] == y,])
  
  
  
  remove_cells <- dplyr::bind_rows(remove_cells)
  
  ind <- apply(remove_cells, 1, function(x) all(is.na(x)))
  remove_cells <- remove_cells[ !ind, ]
  
  remove_cells <- unique(remove_cells[,1])
  
}

match_cols <- function(match_vecs, sv_name){
  out=NULL
  for (i in match_vecs){
    vetor <- i
    vetor <- vetor[vetor != ""]
    key <- data.frame(sample_id=vetor[-1], sv_name=rep(gsub(".*_","", vetor[[1]]), (length(vetor)-1)))  
    out <- rbind(out, key)
  }  
  colnames(out) <- c("sample_id", sv_name) 
  return(out)
}

convert_mt_setting <- function(cell_settings, plot_settings){
  
  test <- readLines(cell_settings)
  
  vecs <- list()
  mtnames <- c()
  for (i in test){
    if (!grepl("#", i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[1]]
      vecs <- append(vecs, lbline)
      mtnames <- append(mtnames, d)    
    }
    
  }
  
  pfx <- unique(gsub("_.*", "", mtnames[grep("_", mtnames)]))
  pfx <- paste0(pfx, "_")
  test <- list()
  for (i in pfx){
    test1 <- list(which(startsWith(mtnames, i)))
    names(test1) = tolower(gsub("_", "", i))
    test <- append(test, test1)
  }
  
  sub_vecs <- list()
  vec_names <- list()
  for (i in test){
    test_vec <- vecs[i]
    sub_vecs <- append(sub_vecs, list(test_vec))
  }
  
  # names(sub_vecs) <- vec_names[1:length(sub_vecs)]
  names(sub_vecs) <- names(test)
  
  sub_vecs <- sub_vecs[unlist(lapply(sub_vecs, length) != 0)]
  
  param_dfs <- purrr::map2(sub_vecs, names(sub_vecs), match_cols)
  
  
  
  if (is.list(param_dfs) & length(param_dfs) != 0) {
    param_dfs <- param_dfs %>%
      Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="sample_id"), .) %>% 
      arrange(sample_id)  
    
    dup_cells <- which(duplicated(param_dfs[,1]))
    if (any(dup_cells)){
      print(paste0("cells ", paste(param_dfs$sample_id[dup_cells], collapse = " "), " found duplicated in cell sets! They will be removed from analysis"))
      param_dfs <- param_dfs[-dup_cells,]
    }
    
    rownames(param_dfs) <- param_dfs[,1]
    
  }
  
  remove_cells <- find_remove_cells(plot_settings, param_dfs)
  
  
  return(list("annotation" = param_dfs, "removed_cells" = remove_cells))
  
}

check_rownames <- function(counts){
  if(grepl("ENST.*", rownames(counts))){
    rownames(counts) <- rownames(counts)
  } else{
    rownames(counts) <- counts[,1]
    counts[,1] <- NULL
  } 
  return(counts)
}

data_append_exp_tag <- function(df, exp_tag){
  names(df) <- paste0(exp_tag, names(df))
  return(df)
}

meta_append_exp_tag <- function(df, exp_tag){
  names(df)[1] <- "sample_id"
  df <- dplyr::mutate(df, sample_id = paste0(exp_tag, sample_id))
  return(df)
}

batch_correct <- function(counts, colData, exp_tags){
  # browser()
  counts <- lapply(counts, check_rownames)
  
  names(counts) <- exp_tags
  names(colData) <- exp_tags
  counts <- purrr::map2(counts, exp_tags, data_append_exp_tag)
  colData <- purrr::map2(colData, exp_tags, meta_append_exp_tag)
  
  # scran mnnCorrect
  counts[[1]] <- counts[[1]][rownames(counts[[1]]) %in% rownames(counts[[2]]),]
  counts[[2]] <- counts[[2]][rownames(counts[[2]]) %in% rownames(counts[[1]]),]
  counts2 <- scran::mnnCorrect(as.matrix(counts[[1]]), as.matrix(counts[[2]]))
  # counts2 <- scran::fastMNN(as.matrix(counts[[1]]), as.matrix(counts[[2]]))

  # joint_expression_matrix <- sce$corrected
  # 
  # joint_tsne <- Rtsne(t(joint_expression_matrix), initial_dims=10, theta=0.75,
  #                       check_duplicates=FALSE, max_iter=200, stop_lying_iter=50, mom_switch_iter=50, perplexity = 2)
  # dataset_labels <- factor(sce$batch)
  # plot(joint_tsne$Y[,1], joint_tsne$Y[,2], pch=c(16,1)[dataset_labels], col=rainbow(length(levels(dataset_labels)))[dataset_labels])
  # 
  return(counts2)
}

make_sce <- function(counts, colData, exp_tags){
  browser()
  counts <- lapply(counts, check_rownames)
  
  names(counts) <- exp_tags
  names(colData) <- exp_tags
  counts <- purrr::map2(counts, exp_tags, data_append_exp_tag)
  colData <- purrr::map2(colData, exp_tags, meta_append_exp_tag)
  
  # scran mnnCorrect
  counts[[1]] <- counts[[1]][rownames(counts[[1]]) %in% rownames(counts[[2]]),]
  counts[[2]] <- counts[[2]][rownames(counts[[2]]) %in% rownames(counts[[1]]),]
  counts2 <- scran::fastMNN(as.matrix(counts[[1]]), as.matrix(counts[[2]]))
  
  cor.exp <- tcrossprod(counts2$rotation[c(1,10),], counts2$corrected)
  dim(cor.exp)
  
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


# load data ------------------------------------------------------

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

# scran::mnnCorrect()

counts2 <- batch_correct(counts, colData, exp_tags)

omat <- do.call(cbind, counts)
mat <- do.call(cbind, counts2$corrected)
colnames(mat) <- NULL
sce <- SingleCellExperiment(list(original=omat, corrected=mat))

# saveRDS(sce, "~/sce_mnn.rds")
colData(sce)$Batch <- rep(names(counts),
                          lapply(counts, ncol))
sce

set.seed(100)
# Using irlba to set up the t-SNE, for speed.
osce <- runPCA(sce, ntop=Inf, method="irlba")
osce <- runTSNE(osce, use_dimred="PCA")
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")

set.seed(100)
csce <- runTSNE(sce, use_dimred="MNN")
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")

multiplot(ot, ct, cols=2)

sce <- make_sce(counts, colData, exp_tags)

saveRDS(sce, outfile)

muraro <- readRDS("~/single_cell_pipeline/bin/muraro.rds")
segerstolpe <- readRDS("~/single_cell_pipeline/bin/segerstolpe.rds")
