
#!/usr/bin/Rscript

#load required libraries and functions
#=====================================
library(optparse)


# shl centroid metadata
# "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_branch_cell_info.csv"

#SHL 20170407
default_expr_mat = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/gene_count_matrix.csv"
default_annotation = "~/single_cell_pipeline/scde_input/shl_0407_w_centroids_cell_info.csv"
default_cell_settings <- "~/single_cell_tools/FACS_0407_2017_SHL_input_files/cell_sets_0407_SHL_20180523.csv"
default_plot_settings <- "~/single_cell_tools/FACS_0407_2017_SHL_input_files/plot_setting_0407_SHL_20180212.csv"
default_pseudotime <- "/home/skevin/tmp/5sh733_shCtrl_PT"
default_out = "/home/skevin"

shl <- mget(ls(pattern = "default"))
save(shl, file = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_phenopath_input.rda")

load( "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_phenopath_input.rda")
list2env(shl, globalenv())


# default input files -----------------------------------------------------

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-e", "--expr_mat"), type="character", default=default_expr_mat,
              help="gene expression input filename [default= %default]", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=default_annotation,
              help="metadata about cells in input file [default= %default]", metavar="character"),
  make_option(c("-c", "--cellset"), type="character", default=default_cell_settings,
              help="tab delimited cell settings file [default= %default]", metavar="character"),
  make_option(c("-p", "--plot_settings"), type="character", default=default_plot_settings,
              help="tab delimited plot settings file [default= %default]", metavar="character"), 
  make_option(c("-t", "--pseudotime"), type="character", default=default_pseudotime,
              help="file that lists cell pseudotimes [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=default_out,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-f", "--feature"), type="character", default=NA, 
              help="the feature (gene or transcript) under analysis", metavar="character")

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
suppressMessages(library(BiocParallel))
suppressMessages(library(tidyverse))
suppressMessages(library(gtools))
suppressMessages(library(cataract))
suppressMessages(library(data.table))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(gplots))
suppressMessages(library(tictoc))
suppressMessages(library(scran))
library(MultiAssayExperiment)
suppressPackageStartupMessages(library(SummarizedExperiment))
edb <- EnsDb.Hsapiens.v86
library(phenopath)


# load required functions -------------------------------------------------

make_sce <- function(cnts, colData, exp_tags){
  # browser()
  check_rownames <- function(cnts){
    if(grepl("ENST.*", rownames(cnts))){
      rownames(cnts) <- rownames(cnts)
    } else{
      rownames(cnts) <- cnts[,1]
      cnts[,1] <- NULL
    } 
    return(cnts)
  }
  
  cnts <- lapply(cnts, check_rownames)
  
  data_append_exp_tag <- function(df, exp_tag){
    names(df) <- paste0(exp_tag, names(df))
    return(df)
  }
  
  meta_append_exp_tag <- function(df, exp_tag){
    names(df)[1] <- "sample_id"
    df <- dplyr::mutate(df, sample_id = paste0(exp_tag, sample_id))
    return(df)
  }
  
  names(cnts) <- exp_tags
  names(colData) <- exp_tags
  cnts <- purrr::map2(cnts, exp_tags, data_append_exp_tag)
  colData <- purrr::map2(colData, exp_tags, meta_append_exp_tag)
  if (length(cnts) > 1){
    if (dim(cnts[[1]]) > dim(cnts[[2]])){
      cnts <- cbind(cnts[[2]], cnts[[1]][match(rownames(cnts[[2]]), rownames(cnts[[1]])),])
    } else if (dim(cnts[[1]] < dim(cnts[[2]]))){
      cnts <- cbind(cnts[[1]], cnts[[2]][match(rownames(cnts[[1]]), rownames(cnts[[2]])),])
    } 
  } else {
    cnts <- cnts[[1]]
  }
  colData <- dplyr::bind_rows(colData)
  colData <- colData[order(colData$sample_id),]
  # rownames(colData) <- colData["sample_id"]
  cnts <- cnts[,order(colnames(cnts))]
  cnts <- as.matrix(cnts)
  
  cnts <- cnts[,(colnames(cnts) %in% colData[,1])]
  
  census_experiment <- SingleCellExperiment(assays=list(counts=cnts), colData=colData, rowData=rownames(cnts))
  return(census_experiment)
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
      dplyr::arrange(sample_id)  
    
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

take_input <- function(prompt = question, interactive = FALSE){
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


# load data ---------------------------------------------------------------

print("loading census matrix")
exprs_mat <- read.csv(opt$expr_mat, row.names = 1)
print("loading cell metadata")
annotation = cataract::safe_read(opt$annotation)

mt_settings <- convert_mt_setting(opt$cellset, opt$plot_settings)

if (!is.null(mt_settings$annotation)) {
  exist_cols <- colnames(mt_settings$annotation)[!colnames(mt_settings$annotation) %in% colnames(annotation)]
  exist_cols <- c("sample_id", exist_cols)
  
  mt_settings$annotation <- dplyr::select(mt_settings$annotation, exist_cols)
  # mt_settings$annotation <- mt_settings$annotation[complete.cases(mt_settings$annotation),]
  oldannotation <- as.data.frame(annotation)
  names(oldannotation) <- tolower(names(oldannotation))
  annotation <- left_join(oldannotation, mt_settings$annotation, by = "sample_id")
  rownames(annotation) <- annotation[,1]
}

annotation<- dplyr::mutate(annotation, cov = dplyr::case_when(
  cluster == "shCtrl" ~ 1,
  cluster == "shRB1" ~ -1
))

# remove cells not present in annotation from expression matrix
exprs_mat <- exprs_mat[,colnames(exprs_mat) %in% annotation[["sample_id"]]]


mae <- readRDS("~/single_cell_tools/tmp/GSE48968-GPL13112.rds")


suppressPackageStartupMessages(library(scater))
cts <- assays(experiments(mae)[["gene"]])[["count_lstpm"]]
tpms <- assays(experiments(mae)[["gene"]])[["TPM"]]
phn <- colData(mae)

sce <- newSCESet(countData = cts, 
                 phenoData = new("AnnotatedDataFrame", data = as.data.frame(phn)))
tpm(sce) <- tpms
exprs(sce) <- log2(tpm(sce) + 1)

is_lps_pam <- grepl("LPS|PAM", sce$description)
sce <- sce[, is_lps_pam]

split <- strsplit(as.character(sce$description), "_", fixed = TRUE)
stimulant <- sapply(split, `[`, 1)
time <- sapply(split, `[`, 2)
sce$stimulant <- stimulant
sce$time <- time

suppressPackageStartupMessages(library(biomaRt))
ensembl_gene_ids <- sapply(strsplit(featureNames(sce), ".", fixed = TRUE), `[`, 1)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
            filters = "ensembl_gene_id",
            values = ensembl_gene_ids,
            mart = mart)

fData(sce)$mgi_symbol <- rep(NA, nrow(sce))

mm2 <- match(bm$ensembl_gene_id, ensembl_gene_ids)
fData(sce)$mgi_symbol[mm2] <- bm$mgi_symbol







## ----example-using-expressionset, eval = FALSE-----------------------------

#  fit <- phenopath(sce, sim$x) # 1
#  fit <- phenopath(sce, "x") # 2
fit <- phenopath(sce, annotation$cov) # 3

## ----initialisation-examples, eval = FALSE---------------------------------
fit <- phenopath(sim$y, sim$x, z_init = 4) # 1, initialise to first principal component
#  fit <- phenopath(sim$y, sim$x, z_init = sim$z) # 2, initialise to true values
#  fit <- phenopath(sim$y, sim$x, z_init = "random") # 3, random initialisation


## ----pcaing, fig.show = 'hold'---------------------------------------------
pca_df <- tbl_df(prcomp(sim$y)$x[,1:2]) %>% 
  mutate(x = factor(sim[['x']]), z = sim[['z']])

ggplot(pca_df, aes(x = PC1, y = PC2, color = x)) +
  geom_point() + scale_colour_brewer(palette = "Set1")

ggplot(pca_df, aes(x = PC1, y = PC2, color = z)) +
  geom_point()

## ----see-results, cache=TRUE-----------------------------------------------
fit <- phenopath(phenopath_in$y, phenopath_in$x, elbo_tol = 1e-6, thin = 40)
print(fit)

## ----plot-elbo-------------------------------------------------------------
plot_elbo(fit)

## ----plot-results, fig.show = 'hold', fig.width = 2.5, fig.height = 2.5----
qplot(sim$z, trajectory(fit)) +
  xlab("True z") + ylab("Phenopath z")
qplot(sim$z, pca_df$PC1) +
  xlab("True z") + ylab("PC1")

## ----print-correlation-----------------------------------------------------
cor(sim$z, trajectory(fit))

## ----beta-df, fig.width = 6, fig.height = 3--------------------------------
gene_names <- paste0("gene", seq_len(ncol(fit$m_beta)))
df_beta <- data_frame(beta = interaction_effects(fit),
                      beta_sd = interaction_sds(fit),
                      is_sig = significant_interactions(fit),
                      gene = gene_names)

df_beta$gene <- fct_relevel(df_beta$gene, gene_names)

ggplot(df_beta, aes(x = gene, y = beta, color = is_sig)) + 
  geom_point() +
  geom_errorbar(aes(ymin = beta - 2 * beta_sd, ymax = beta + 2 * beta_sd)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x = element_blank()) +
  ylab(expression(beta)) +
  scale_color_brewer(palette = "Set2", name = "Significant")

## ----graph-largest-effect-size---------------------------------------------
which_largest <- which.max(df_beta$beta)

df_large <- data_frame(
  y = sim[['y']][, which_largest],
  x = factor(sim[['x']]),
  z = sim[['z']]
)

ggplot(df_large, aes(x = z, y = y, color = x)) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  stat_smooth()

## ----construct-sceset, warning = FALSE-------------------------------------

## ----cavi-tuning, eval = FALSE---------------------------------------------
#  fit <- phenopath(sim$y, sim$x,
#                   maxiter = 1000, # 1000 iterations max
#                   elbo_tol = 1e-2, # consider model converged when change in ELBO < 0.02%
#                   thin = 20 # calculate ELBO every 20 iterations
#                   )

## ----sessioninfo-----------------------------------------------------------
sessionInfo()
