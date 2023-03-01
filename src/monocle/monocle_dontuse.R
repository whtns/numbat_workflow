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
default_cell_info <- "~/single_cell_tools/FACS_0407_2017_SHL_input_files/cell_sets_0407_SHL_20180523.csv"
default_plot_settings <- "~/single_cell_tools/FACS_0407_2017_SHL_input_files/plot_setting_0407_SHL_20180212.csv"
default_out <- "~/1stdata"

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
  make_option(c("-c", "--cellset"), type="character", default=default_cell_info,
              help="tab delimited cell settings file [default= %default]", metavar="character"),
  make_option(c("-p", "--plot_settings"), type="character", default=default_plot_settings,
              help="tab delimited plot settings file [default= %default]", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NA,
              help="output file name [default= %default]", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  missing_args <- paste(names(which(sapply(opt, is.na))), collapse = " ")
  stop(paste0("Please provide all necessary arguments. Missing: ", missing_args), call.=FALSE)
}

# load required libraries -------------------------------------------------
if (!grepl(".pdf", opt$outfile)){
  print("output file does not have .pdf extension, appending...")
  out_pdf <- path.expand(paste0(opt$outfile, ".pdf"))
} else {
  out_pdf <- path.expand(opt$outfile)
}


suppressMessages(library(monocle))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(SCnorm))
suppressMessages(library(tidyverse))
suppressMessages(library(gtools))
suppressMessages(library(tictoc))
suppressMessages(library(memoise))
suppressMessages(library(Biobase))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

fc <- cache_filesystem("~/.rcache")



# load required functions -------------------------------------------------

take_input <- function(prompt = question, interactive = TRUE){
  if(interactive){
    param <- readline(prompt=question)
  } else {
    param <- readLines("stdin", n=1)
  }
  return(param)
}

lookup_genes <- function(txname){
  
  txs <- transcripts(edb, filter = TxIdFilter(txname),
                     columns = c("symbol"))
  return(txs$symbol)
  
}

lookup_transcripts <- function(tx_df){
  
  txids <- rownames(tx_df)
  
  txs <- transcripts(edb, filter = TxIdFilter(txids), columns = c("symbol"))
  
  names(txids) <- txs$symbol
  
  return(txids)
}



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

batch_plot_genes_in_pseudotime <- function(colorval, ordering_genes, HSMM){
  for (i in colorval){
    plot_vec = c(seq(1, length(ordering_genes), 5), length(ordering_genes))
    for (m in seq(plot_vec[-1])){
      beg_plot = plot_vec[m]
      end_plot = plot_vec[m+1]
      nplots = end_plot - beg_plot + 1
      print(plot_genes_in_pseudotime(HSMM[ordering_genes[plot_vec[m]:plot_vec[m+1]],],
                                     color_by = i, label_by_short_name = FALSE))
    }
  }
}

batch_plot_multiple_branches_pseudotime <- function(colorval, top_BEAM, top_genes, branches, HSMM){
  
  for (i in colorval){
    plot_vec = c(seq(1, top_genes, 5), top_genes)
    for (m in seq(plot_vec[-1])){
      beg_plot = plot_vec[m]
      end_plot = plot_vec[m+1]
      nplots = end_plot - beg_plot + 1
      if (nplots > 3){
        ncols <- 3
      } else {
        ncols = nplots
      }
      print(plot_multiple_branches_pseudotime(HSMM[rownames(top_BEAM)[plot_vec[m]:plot_vec[m+1]],],
                                              branches=branches,
                                              color_by = 'Branch',
                                              label_by_short_name = FALSE,
                                              nrow = ceiling(nplots/3), ncol = ncols))
      
    }
    
  }
}

batch_plot_genes_branched_pseudotime <- function(colorval, top_BEAM, top_genes, branch_states, HSMM){
  
  for (i in colorval){
    plot_vec = c(seq(1, top_genes, 5), top_genes)
    for (m in seq(plot_vec[-1])){
      beg_plot = plot_vec[m]
      end_plot = plot_vec[m+1]
      nplots = end_plot - beg_plot + 1
      if (nplots > 3){
        ncols <- 3
      } else {
        ncols = nplots
      }
      print(plot_genes_branched_pseudotime(HSMM[rownames(top_BEAM)[plot_vec[m]:plot_vec[m+1]],],
                                           branch_states = branch_states, 
                                           color_by = 'Branch', 
                                           label_by_short_name = FALSE,
                                           method = 'loess', 
                                           nrow = ceiling(nplots/3), ncol = ncols))
    }
    
  }
}

GM_state <- function(cds, pt_param = "day", pt_paramval = "day_0", branches = NULL){
  
  if (!is.null(branches)){
    cds <- cds[,pData(cds)$State %in% branches]
  }
  
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)[,pt_param])[,pt_paramval]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}




run_qc <- function(sce, param, paramval, colorval, annotation=NULL){
  
  
  # out_pdf <- path.expand(out_pdf)
  # out_pdf <- gsub(".pdf", "_census.pdf", out_pdf)
  
  # pdf(out_pdf)
  # subset summarized experiment according to user conditions
  
  if(grepl(",", paramval)){
    paramval <- unlist(strsplit(paramval, ","))
  }
  
  if(param == '' | paramval == ''){
    subset_sce <- sce
  } else{
    subset_sce <- sce[, (sce[[param]] %in% paramval)]
  }
  
  if(grepl(",", colorval)){
    colorval <- unlist(strsplit(colorval, ","))
  }
  
  sample_sheet <- as.data.frame(SummarizedExperiment::colData(subset_sce))
  gene_short_names <- make.unique(lookup_genes(SummarizedExperiment::rowData(subset_sce)[,1]))
  transcripts_annotation <- data.frame("transcripts" = SummarizedExperiment::rowData(subset_sce), "gene_short_name" = gene_short_names)
  rownames(transcripts_annotation) <- transcripts_annotation[,1]
  
  pd <- new("AnnotatedDataFrame", data=sample_sheet)
  fd <- new("AnnotatedDataFrame", data=transcripts_annotation)
  
  # First create a CellDataSet from the relative expression levels
  HSMM <- newCellDataSet(as.matrix(assays(subset_sce)$counts),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=0.1,
                         expressionFamily=tobit(Lower=0.1))
  # Next, use it to estimate RNA counts
  rpc_matrix <- relative2abs(HSMM)
  
  NA_count =sum(is.na(rpc_matrix))
  rpc_matrix <- na.replace(rpc_matrix, 0.)
  
  OnlyZeros = (colSums(rpc_matrix)==0)
  paste(sum(OnlyZeros), "cells have zero reads in total, and there were", NA_count, "NA values before replacement to NA -> 0")
  
  Valid = which(!OnlyZeros)
  rpc_matrix = rpc_matrix[ , Valid ]; dim(rpc_matrix)
  
  if (length(which(OnlyZeros)) >0){
    pd <- pd[-which(OnlyZeros),]  
  }
  
  
  # Now, make a new CellDataSet using the RNA counts
  
  HSMM <- newCellDataSet(as(rpc_matrix, "sparseMatrix"),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=negbinomial.size())
  
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  # filter by gene expression
  # # keep genes that are expressed in at least 5 cells
  # min_expression <- 0.1
  # HSMM <- detectGenes(HSMM, min_expr=min_expression)
  # expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 5))
  
  # filter by gene expression (default monocle settings)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  expressed_genes <- row.names(subset(fData(HSMM),
                                      num_cells_expressed >= 10))
  
  # look at distribution of mRNA totals across cells
  pData(HSMM)$Total_mRNAs <- Matrix::colSums(Biobase::exprs(HSMM))
  
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
  
  upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                       2*sd(log10(pData(HSMM)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                       2*sd(log10(pData(HSMM)$Total_mRNAs)))
  
  # remove cells outside safe range of plot ---------------------------------
  
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
                 pData(HSMM)$Total_mRNAs < upper_bound]
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  # verify lognormal distribution of expression values ----------------------
  
  # Log-transform each value in the expression matrix.
  L <- log(Biobase::exprs(HSMM[expressed_genes,]))
  
  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  
  # Plot the distribution of the standardized gene expression values.
  exp_dist <- qplot(value, geom = "density", data = melted_dens_df) +
    stat_function(fun = dnorm, size = 0.5, color = 'red') +
    xlab("Standardized log(FPKM)") +
    ylab("Density")
  print(exp_dist)
  
  # cluster cells without marker genes --------------------------------------
  
  # develop list of genes used for clustering
  disp_table <- dispersionTable(HSMM)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  print(plot_ordering_genes(HSMM))
  
  # plot variance explained by each principle component
  # HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
  print(plot_pc_variance_explained(HSMM, return_all = F)) # norm_method='log'
  
  # cluster cells by tsne
  HSMM_red <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                              reduction_method = 'tSNE', verbose = T)
  HSMM_red <- clusterCells(HSMM_red, num_clusters = 2)
  
  for (i in colorval){
    print(plot_cell_clusters(HSMM_red, 1, 2, color = i))    
  }
  
  
  # subtract uninteresting sources of variation (batch, etc.)
  # HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
  #                         reduction_method = 'tSNE',
  #                         residualModelFormulaStr = "~index_i7 + num_genes_expressed",
  #                         verbose = T)
  # HSMM <- clusterCells(HSMM, num_clusters = 2)
  # print(plot_cell_clusters(HSMM, 1, 2, color = "treatment_group")
  
  return(list("HSMM" = HSMM, "expressed_genes" = expressed_genes))  
  
}

run_census <- function(sce, param, paramval, colorval, annotation=NULL){
  
  # out_pdf <- path.expand(out_pdf)
  # out_pdf <- gsub(".pdf", "_census.pdf", out_pdf)
  
  # pdf(out_pdf)
  # subset summarized experiment according to user conditions
  
  if(grepl(",", paramval)){
    paramval <- unlist(strsplit(paramval, ","))
  }
  
  if(param == '' | paramval == ''){
    subset_sce <- sce
  } else{
    subset_sce <- sce[, (sce[[param]] %in% paramval)]
  }
  
  if(grepl(",", colorval)){
    colorval <- unlist(strsplit(colorval, ","))
  }
  
  sample_sheet <- as.data.frame(SummarizedExperiment::colData(subset_sce))
  gene_short_names <- make.unique(lookup_genes(SummarizedExperiment::rowData(subset_sce)[,1]))
  transcripts_annotation <- data.frame("transcripts" = SummarizedExperiment::rowData(subset_sce), "gene_short_name" = gene_short_names)
  rownames(transcripts_annotation) <- transcripts_annotation[,1]
  
  pd <- new("AnnotatedDataFrame", data=sample_sheet)
  fd <- new("AnnotatedDataFrame", data=transcripts_annotation)
  
  # First create a CellDataSet from the relative expression levels
  HSMM <- newCellDataSet(as.matrix(assays(subset_sce)$counts),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=0.1,
                         expressionFamily=tobit(Lower=0.1))
  # Next, use it to estimate RNA counts
  rpc_matrix <- relative2abs(HSMM)
  
  NA_count =sum(is.na(rpc_matrix))
  rpc_matrix <- na.replace(rpc_matrix, 0.)
  
  OnlyZeros = (colSums(rpc_matrix)==0)
  paste(sum(OnlyZeros), "cells have zero reads in total, and there were", NA_count, "NA values before replacement to NA -> 0")
  
  Valid = which(!OnlyZeros)
  rpc_matrix = rpc_matrix[ , Valid ]; dim(rpc_matrix)
  
  if (length(which(OnlyZeros)) >0){
    pd <- pd[-which(OnlyZeros),]  
  }
  
  
  # Now, make a new CellDataSet using the RNA counts
  
  HSMM <- newCellDataSet(as(rpc_matrix, "sparseMatrix"),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=negbinomial.size())
  
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  # filter by gene expression
  # # keep genes that are expressed in at least 5 cells
  # min_expression <- 0.1
  # HSMM <- detectGenes(HSMM, min_expr=min_expression)
  # expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 5))
  
  # filter by gene expression (default monocle settings)
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  expressed_genes <- row.names(subset(fData(HSMM),
                                      num_cells_expressed >= 10))
  
  # look at distribution of mRNA totals across cells
  pData(HSMM)$Total_mRNAs <- Matrix::colSums(Biobase::exprs(HSMM))
  
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
  
  upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                       2*sd(log10(pData(HSMM)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                       2*sd(log10(pData(HSMM)$Total_mRNAs)))
  
  # remove cells outside safe range of plot ---------------------------------
  
  HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
                 pData(HSMM)$Total_mRNAs < upper_bound]
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  # verify lognormal distribution of expression values ----------------------
  
  # Log-transform each value in the expression matrix.
  L <- log(Biobase::exprs(HSMM[expressed_genes,]))
  
  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  
  
  # cluster cells without marker genes --------------------------------------
  
  # develop list of genes used for clustering
  disp_table <- dispersionTable(HSMM)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  
  # cluster cells by tsne
  HSMM_red <- reduceDimension(HSMM, max_components = 2, num_dim = 6,
                              reduction_method = 'tSNE', verbose = T)
  HSMM_red <- clusterCells(HSMM_red, num_clusters = 2)
  
  # subtract uninteresting sources of variation (batch, etc.)
  # HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 2,
  #                         reduction_method = 'tSNE',
  #                         residualModelFormulaStr = "~index_i7 + num_genes_expressed",
  #                         verbose = T)
  # HSMM <- clusterCells(HSMM, num_clusters = 2)
  # print(plot_cell_clusters(HSMM, 1, 2, color = "treatment_group")
  
  return(list("HSMM" = HSMM, "expressed_genes" = expressed_genes))  
  
}
mrun_census <- memoise::memoise(run_census, cache = fc)

run_monocle <- function(census_output, compareval, reduceval, colorval, start_cells = FALSE, pt_param = "day", pt_paramval = "day_0", diff_by_var = FALSE, gene_list = NULL){
   if(grepl(",", colorval)){
    colorval <- unlist(strsplit(colorval, ","))
  }
  
  # check if gene list file can be found

  HSMM <- census_output[["HSMM"]]
  expressed_genes <- census_output[["expressed_genes"]]
  
  # subset summarized experiment according to user conditions
  
  # choose genes that define a cell's progress ------------------------------
  
  # use cell metadata -------------------------------------------------------
  
  if (diff_by_var){
    
    print("running differential expression test")
    tic("finished differentiial expression with")
    diff_test_res <- differentialGeneTest(HSMM[expressed_genes,], fullModelFormulaStr = paste0("~", compareval),
    reducedModelFormulaStr = paste0("~", reduceval), cores = 6)
    # 
    # diff_test_res <- readRDS("~/diff_test_res.rds")
    
    toc()   
    
    if(!is.null(gene_list)){
      # browser()
      gene_list <- read.csv(gene_list, header = FALSE)
      ordering_genes <- rownames(diff_test_res)[which(gene_list[,1] %in% rownames(diff_test_res))][1:1000]
    } else {
      ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
      
      #select top 1000 signif. genes
      # ordering_genes <-  row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
    }
    
    
    # set differentially expressed genes in the HSMM object
    HSMM <- setOrderingFilter(HSMM, ordering_genes)
    
    
  } else { #unsupervised differential expression
    
    # use dpFeature -----------------------------------------------------------
    
    
    HSMM <- detectGenes(HSMM, min_expr = 0.1)
    fData(HSMM)$use_for_ordering <- fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)    
    
    print(plot_pc_variance_explained(HSMM, return_all = F))
    
    HSMM_red <- reduceDimension(HSMM, max_components = 2,
                                norm_method = 'log',
                                num_dim = 3,
                                reduction_method = 'tSNE',
                                verbose = T)
    
    HSMM_red <- clusterCells(HSMM_red, verbose = F)
    
    # check clustering results
    print(plot_cell_clusters(HSMM_red, color_by = 'as.factor(Cluster)'))
    
    for (i in colorval){
      print(plot_cell_clusters(HSMM_red, color_by = paste0('as.factor(', i, ')')))      
    }
    
    
    # provide decision plot
    print(plot_rho_delta(HSMM_red, rho_threshold = 2, delta_threshold = 4 ))
    
    # rerun based on user-defined threshold
    HSMM_red <- clusterCells(HSMM_red, rho_threshold = 2,
                             delta_threshold = 4,
                             skip_rho_sigma = T,
                             verbose = F)
    
    # check final clustering
    print(plot_cell_clusters(HSMM_red, color_by = 'as.factor(Cluster)'))
    
    for (i in colorval){
      print(plot_cell_clusters(HSMM_red, color_by = paste0('as.factor(', i, ')')))      
    }
    
    
    
    # perform differential expression -----------------------------------------
    
    # find expressed genes    
    HSMM_expressed_genes <-  row.names(subset(fData(HSMM), num_cells_expressed >= 10))
    
    print("running differential expression test")
    tic("finished differentiial expression with")
    diff_test_res <- differentialGeneTest(HSMM_red[HSMM_expressed_genes,],
                                                 fullModelFormulaStr = '~Cluster',
                                                 cores = 6)
    toc()
    
    if(!is.null(gene_list)){
      gene_list <- read.csv(gene_list, header = FALSE)
      ordering_genes <- rownames(diff_test_res)[which(gene_list[,1] %in% rownames(diff_test_res))][1:1000]
    } else {
      #select top 1000 signif. genes
      ordering_genes <-  row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
      
    }
    
    
    
    
    HSMM <- setOrderingFilter(HSMM, ordering_genes = ordering_genes)
  }
  
  HSMM <-   reduceDimension(HSMM, method = 'DDRTree')
  
  HSMM <-  orderCells(HSMM)
  
  for (i in colorval){
    print(plot_cell_trajectory(HSMM, color_by = i, show_cell_names=FALSE, cell_size = 4, show_tree = TRUE, show_backbone = FALSE) + geom_text(aes(label = sample_name), size = 1))
    
    print(plot_cell_trajectory(HSMM, color_by = i) +
            facet_wrap(~State, nrow = 1))    
  }
  
  
  # plot trajectory colored by pseudotime -----------------------------------
  if(start_cells){
    
    HSMM <- orderCells(HSMM, root_state =  GM_state(HSMM, pt_param, pt_paramval))
    print(plot_cell_trajectory(HSMM, color_by = "Pseudotime", show_cell_names = TRUE, cell_name_size = 1, cell_size =4, show_tree = TRUE))
    
    print(plot_cell_trajectory(HSMM, color_by = "Pseudotime", show_cell_names = TRUE, cell_name_size = 1, cell_size =4, show_tree = TRUE) +
            facet_wrap(~State, nrow = 1)) 
  }
  
  png(file = "test.png")
  
  pt_hmap <- plot_pseudotime_heatmap(HSMM[ordering_genes,],
                                     num_clusters = 3,
                                     cores = 6,
                                     show_rownames = T,
                                     use_gene_short_name = F, return_heatmap = TRUE)
  dev.off()
  file.remove('test.png')
  # make font size of row annottations smaller
  
  pt_hmap$gtable$grobs[[3]]$gp$fontsize <- 2
  
  # browser()
  return_trx <- pt_hmap$gtable$grobs[[3]]$label
  return_genes <- lookup_genes(return_trx)
  hmap_genes <- data.frame("ensembl_transcript_id" = return_trx, "gene_symbol" = return_genes, "pval" = diff_test_res[return_trx,]$pval, "qval" = diff_test_res[return_trx,]$qval)
  hmap_genes_csv <- gsub(".pdf", paste0("_beam_heatmap_genes.csv"), out_pdf)
  write.csv(hmap_genes, hmap_genes_csv)
  
  hmap_genes <- dplyr::mutate(hmap_genes, new_labels = paste(ensembl_transcript_id, gene_symbol, sep = "\t"))
  
  pt_hmap$gtable$grobs[[3]]$label <- hmap_genes$new_labels
  
  gridExtra::grid.arrange(pt_hmap[[4]])
  
  
  # batch_plot_genes_in_pseudotime(colorval, ordering_genes, HSMM)
  print(paste0("saving monocle output: ", monocle_pdf))
  return(list("HSMM" = HSMM, "diff_test_res" = diff_test_res))
}



mrun_monocle <- memoise::memoise(run_monocle, cache = fc)

run_BEAM <- function(HSMM, branch_point, branches, pt_param, pt_paramval, colorval, top_genes){
  
  if(grepl(",", colorval)){
    colorval <- unlist(strsplit(colorval, ","))
  }
  
  # check if branches present in HSMM
  
  root_branch <- GM_state(HSMM, pt_param, pt_paramval, branches)
  branch_states = branches[which(branches != root_branch)]
  
  # test <- buildBranchCellDataSet(HSMM, progenitor_method = "sequential_split", branch_point = 2,
  #                                branch_labels = NULL, stretch = TRUE)
  # 
  # test <- buildBranchCellDataSet(HSMM, progenitor_method = "duplicate", branch_point = 2,
  #                                branch_labels = NULL, stretch = TRUE)
  
  BEAM_res <- BEAM(HSMM, branch_point = branch_point, cores = 6,  progenitor_method = "sequential_split")
  
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("pval", "qval")]
  
  return(BEAM_res)

}

run_hmap <- function(HSMM, BEAM_res, branches, pt_param, pt_paramval, colorval, top_genes, out_pdf){
  beam_pdf <- gsub(".pdf", "_beam.pdf", out_pdf)
  pdf(beam_pdf, height = 10)
  # browser()
  plotted_trx <- row.names(subset(BEAM_res, qval < 1e-4))
  gene_df <- data.frame(row.names = plotted_trx, gene_symbol = lookup_genes(plotted_trx))
  
  png(file = "test.png")
  mult_hmap <- plot_multiple_branches_heatmap(HSMM[plotted_trx,],
                                              branches = branches,
                                              num_clusters = 2,
                                              cores = 6,
                                              use_gene_short_name = F,
                                              show_rownames = T, return_heatmap = TRUE)
  dev.off()
  
  mult_hmap$gtable$grobs[[3]]$gp$fontsize <- 4
  
  return_trx <- mult_hmap$gtable$grobs[[3]]$label
  return_genes <- lookup_genes(return_trx)
  
  hmap_genes <- data.frame("ensembl_transcript_id" = return_trx, "gene_symbol" = return_genes, "pval" = BEAM_res[return_trx,]$pval, "qval" = BEAM_res[return_trx,]$qval)
  hmap_genes_csv <- gsub(".pdf", paste0("_beam_heatmap_genes.csv"), out_pdf)
  write.csv(hmap_genes, hmap_genes_csv)
  
  hmap_genes <- dplyr::mutate(hmap_genes, new_labels = paste(ensembl_transcript_id, gene_symbol, sep = "\t"))
  
  mult_hmap$gtable$grobs[[3]]$label <- hmap_genes$new_labels
  
  gridExtra::grid.arrange(mult_hmap[[4]])
  
  
  root_branch <- GM_state(HSMM, pt_param, pt_paramval)
  branch_states = branches[-root_branch]
  
  # plot individual genes ---------------------------------------------------
  # top_BEAM <- BEAM_res[1:top_genes,]
  top_BEAM <- BEAM_res[return_trx,]
  top_genes <- length(return_trx)
  
  
  
  batch_plot_multiple_branches_pseudotime(colorval, top_BEAM, top_genes, branches, HSMM)
  # 
  # batch_plot_genes_branched_pseudotime(colorval, top_BEAM, top_genes, branch_states, HSMM)
  
  dev.off()
  return(mult_hmap)
  
}


# load data ---------------------------------------------------------------

# load summarized experiment
print("loading data")

sce <- readRDS(opt$infile)

mt_settings <- convert_mt_setting(opt$cellset, opt$plot_settings)

annotation <- mt_settings$annotation

if (!is.null(annotation)) {
  exist_cols <- colnames(annotation)[!colnames(annotation) %in% colnames(colData(sce))]
  exist_cols <- c("sample_id", exist_cols)
  
  annotation <- dplyr::select(annotation, exist_cols)
  # annotation <- annotation[complete.cases(annotation),]
  oldcoldata <- as.data.frame(colData(sce))
  names(oldcoldata) <- tolower(names(oldcoldata))
  newcoldata <- left_join(oldcoldata, annotation, by = "sample_id")
  
  sce <- SingleCellExperiment(assays=list(counts=counts(sce)), colData=newcoldata, rowData=rownames(counts(sce)))
  sce <- sce[, !sce$sample_id %in% mt_settings$removed_cells]
}


while (TRUE) {
  # if (ctr < 3){
  #   print(question)
  # }
  question = "Choose from following:\n[S] Show Cell Info\n[Q] Run QC\n[M] Run Monocle\n[B] Run BEAM (differential expression)\n[X] Exit "  
  cat(question)
  # action <- toupper(readLines("stdin",n=1))
  action <- readline(prompt=question)
  
  if(action == "X"){
    break
  } else if(action == "S"){
    print("showing cell info: ")
    print(colData(sce))
  } else if(action == "Q"){
    
    # run QC ---------------------------------
    question <- "Enter parameter to subset on (day, treatment_group, cluster, etc.); blank if none: "
    cat(question)
    param <- take_input(prompt=question)
    
    question <- paste0("Enter subset values; comma separated if multiple; blank if none (", paste(levels(sce[[param]]), collapse = ','), "): ")
    cat(question)
    paramval <- take_input(prompt=question)
    
    question <- "Enter parameter on which to color (ex. day, treatment_group, etc.); comma separated if multiple "
    cat(question)
    colorval <- take_input(prompt=question)
    
    census_pdf <- gsub(".pdf", "_census.pdf", out_pdf)
    
    HSMM_rds <- gsub(".pdf", "_hsmm.rds", out_pdf)
    
    pdf(census_pdf)
    census_output <- run_qc(sce, param, paramval, colorval)  
    dev.off()
    print(paste0("saving census output: ", census_pdf))
    
    
  } else if (action == "M"){

    forget(mrun_monocle)
    
    question <- "use prior data? (y,n) "
    cat(question)
    use_prior <- take_input(prompt=question)
    use_prior  <- ifelse(use_prior == "y", TRUE, FALSE)
    
    # run Monocle ---------------------------------
    question <- "Enter parameter to subset on (day, treatment_group, cluster, etc.); blank if none: "
    cat(question)
    param <- take_input(prompt=question)
    
    question <- paste0("Enter subset values; comma separated if multiple; blank if none (", paste(levels(sce[[param]]), collapse = ','), "): ")
    cat(question)
    paramval <- take_input(prompt=question)
    
    question <- "Enter parameter on which to color (ex. day, treatment_group, etc.); comma separated if multiple "
    cat(question)
    colorval <- take_input(prompt=question)
    
    question <- "start cells known? (y,n) "
    cat(question)
    start_cells <- take_input(prompt=question)
    start_cells  <- ifelse(start_cells == "y", TRUE, FALSE)
    
    if (start_cells){
      question <- "set pseudotime zero with what param? "
      cat(question)
      pt_param <- take_input(prompt=question)
    } else {
      pt_param <- "day"
    }
    
    if (start_cells){
      question <- question <- "set pseudotime zero with what value? "
      cat(question)
      pt_paramval <- take_input(prompt=question)
    } else {
      pt_paramval <- "day_0"
    }
    
    question <- "supply list of genes for differential expression (as .csv or .txt) "
    cat(question)
    gene_list <- take_input(prompt=question)
    if(gene_list == ""){
      gene_list = NULL
    } else {
      gene_list_pat <- gene_list
      gene_list <- tryCatch(list.files("~/single_cell_tools/", pattern = gene_list, full.names = TRUE, recursive = TRUE)[[1]], 
               error=function(e) stop(paste0("file ", gene_list_pat, " can't be found")))
      print(paste0("found file listing genes correlated with pcs: ", gene_list))
      
    }

    question <- "use metadata to inform trajectory? (y,n) "
    cat(question)
    meta_use <- take_input(prompt=question)
    meta_use  <- ifelse(meta_use == "y", TRUE, FALSE)
    
    if (meta_use){
      question <- "Enter parameter on which to run differential expression (day, treatment_group, etc.) "
      cat(question)
      compareval <- take_input(prompt=question)
    } else { 
      compareval <- ""
    }
    
    if (meta_use){
      question <- "Enter parameter to compare for Likelihood Ratio Test (LRT) "
      cat(question)
      reduceval <- take_input(prompt=question)
      if (reduceval == ''){
        reduceval <- "1"
      }
    } else { 
      reduceval <- "1"
    }
    if(!use_prior){
      forget(mrun_census)
      census_args <- list(sce, param, paramval, colorval)
      saveRDS(census_args, file = "~/.census_args.rds")
    } else {
      census_args <- readRDS("~/.census_args")
    }
    
    census_output <- do.call(mrun_census, census_args) 
    
    monocle_args <- list(census_output, compareval, reduceval, colorval, start_cells, pt_param, pt_paramval, diff_by_var = meta_use, gene_list)
    
    saveRDS(monocle_args, file = "~/.monocle_args.rds")
    
    monocle_pdf <- gsub(".pdf", "_monocle.pdf", out_pdf)
    pdf(monocle_pdf)
    
    monocle_output <- do.call(mrun_monocle, monocle_args)
    
    dev.off()
    HSMM <- monocle_output$HSMM
    diff_test_res <- monocle_output$diff_test_res
    
    
  } else if (action == "B"){
    
    question <- "use prior data? (y,n) "
    cat(question)
    use_prior <- take_input(prompt=question)
    use_prior  <- ifelse(use_prior == "y", TRUE, FALSE)
    
    if(!use_prior){
      forget(mrun_census)
      forget(mrun_monocle)
    } else {
      
      census_args <- readRDS("~/.census_args.rds")
      monocle_args <- readRDS("~/.monocle_args.rds")
      
      census_output <- do.call(mrun_census, census_args) 
      monocle_output <- do.call(mrun_monocle, monocle_args)
      HSMM <- monocle_output$HSMM
      diff_test_res <- monocle_output$diff_test_res
    }
    
    question <- "which branchpoint? "
    cat(question)
    branch_point <- take_input(prompt=question)
    branch_point <- as.numeric(branch_point)
    
    question <- "which branches do you want to compare (1,2,3; etc.); split values by , "
    cat(question)
    branches <- take_input(prompt=question)
    branches <- as.numeric(unlist(strsplit(branches, ",")))
    
    question <- "Enter parameter on which to color (ex. day, treatment_group, etc.); comma separated if multiple "
    cat(question)
    colorval <- take_input(prompt=question)
    
    # question <- "number of genes to check correlation "
    # cat(question)
    # top_genes <- as.numeric(take_input(prompt=question))
    # # top_genes <- read.csv(top_genes, header = FALSE)[,1]
    
    if (identical(HSMM@auxOrderingData$DDRTree$branch_points, character(0))){
      print("no branches in supplied pseudotime")
    } else {
      
      
      beam_pdf <- gsub(".pdf", "_beam.pdf", out_pdf)
      print("running BEAM differential expression ")
      tic("beam completed with")
      BEAM_res <- run_BEAM(HSMM, branch_point, branches, pt_param, pt_paramval, colorval)
      pdf(beam_pdf, height = 10)
      mult_hmap <- run_hmap(HSMM, BEAM_res, branches, pt_param, pt_paramval, colorval, top_genes = 3, out_pdf)
      
      
      
      dev.off()
      toc()
    }
    
  } else if (action == "I"){
    browser()
    print("a")
  } 
}

