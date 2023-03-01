#!/usr/bin/Rscript

#load required libraries and functions
#=====================================
library(optparse)


# shl centroid metadata
# "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_branch_cell_info.csv"

# shl <- mget(ls(pattern = "default"))
# save(shl, file = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_0407_deseq_input.rda")

load( "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_0407_deseq_input.rda")
list2env(shl, globalenv())


# default input files -----------------------------------------------------

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=default_infile,
              help="gene expression input filename [default= %default]", metavar="character"),
  make_option(c("-c", "--cell_info"), type="character", default=default_cell_info,
              help="metadata about cells in input file [default= %default]", metavar="character"),
  make_option(c("-g", "--groupfile"), type="character", default=default_groupfile,
              help="file that defines comparison groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=default_out,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-f", "--feature"), type="character", default="transcript", 
              help="the feature (gene or transcript) under analysis", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  missing_args <- paste(names(which(sapply(opt, is.na))), collapse = " ")
  stop(paste0("Please provide all necessary arguments. Missing: ", missing_args), call.=FALSE)
}

run_edger <- function(summarized_exp, full_model){
  full_model <- gsub("~", "", full_model)
  dgelist <- DGEList(counts=summarized_exp$counts, group=summarized_exp$sample_sheet[[full_model]])
  
  design <- model.matrix(~0+group, data=dgelist$samples)
  colnames(design) <- levels(dgelist$samples$group)
  
  dgelist <- estimateDisp(dgelist, design, robust=TRUE)
  
  fit <- glmQLFit(dgelist, design, robust=TRUE)
  
  qlf <- glmQLFTest(fit, coef=1:ncol(design))
  
  top <- topTags(qlf, n=10000L, p.value = 0.01)
  
  test <- edgeR::cpm(dgelist)[rownames(top),]
  
  baseMeanPerLvl <- sapply( levels(dgelist$samples$group), function(lvl) rowMeans( test[,rownames(dgelist$samples[dgelist$samples$group == lvl,])] ) )
  colnames(baseMeanPerLvl) <- paste0("CPM_", colnames(baseMeanPerLvl))
  
  num_groups <- length(levels(dgelist$samples$group))
  
  top$table[,1:num_groups] <- baseMeanPerLvl[,1:num_groups]
  colnames(top$table)[1:num_groups] <- colnames(baseMeanPerLvl)[1:num_groups] 
  
  # factor(mat%*%(1:5), labels = colnames(mat)
  
  return(top)
  
}


# load required libraries -------------------------------------------------
print("loading required libraries")
suppressMessages(library(BiocParallel))
suppressMessages(library(DESeq2))
suppressMessages(library(tidyverse))
suppressMessages(library(gtools))
suppressMessages(library(cataract))
suppressMessages(library(data.table))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(gplots))
suppressMessages(library(shiny))
suppressMessages(library(heatmaply))
suppressMessages(library(shinyHeatmaply))
edb <- EnsDb.Hsapiens.v86
suppressMessages(library(edgeR))


#Generate RNA-seq matrix
#Set parameters
#=================================================

# make output direcotry ---------------------------------------------------
suppressWarnings(dir.create(opt$out))

user_groups <- read.table(textConnection(gsub(" ", "\t", readLines(opt$groupfile))), sep="\t", header = TRUE, stringsAsFactors = FALSE, fill=TRUE)

counts <- cataract::safe_read(opt$infile)

sample_sheet <- read.csv(opt$cell_info, row.names=1) 
# 
# trace("prep_counts_and_sample_sheet", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))

summarized_exp <- purrr::map2(user_groups$subset_parameter, user_groups$subset_param_values, cataract::prep_counts_and_sample_sheet, counts, sample_sheet, strip_col = TRUE)

# continue with deseq2 analysis -------------------------------------------

# full_list <- list()
# for (row in 1:nrow(user_groups)){
#     group_list <- as.list(user_groups[row,])
#     names(group_list) <- names(user_groups)
#     print(row)
#     full_list[[row]] <- group_list
#   }  

full_model_params <- gsub("[[:punct:]]", "", user_groups$full_model)
reduced_model_params <- gsub(":[[punct]]:", "", user_groups$reduced_model) 


# factorize all model parameters ------------------------------------------
summarized_exp <- purrr::map2(summarized_exp, full_model_params, function(x,y){ x[["sample_sheet"]][[y]] <- as.factor(x[["sample_sheet"]][[y]]); return(x)})

#specify model design
#=============================

question = "Choose from following:\n[L] Run Likelihood Ratio Test\n[D] Run Group Comparison\n[X] Exit\n "

while (TRUE) {
  # if (ctr < 3){
  #   print(question)
  # }
  
  cat(question)
  action <- toupper(readLines("stdin",n=1))
  # action <- readline(prompt=question)
  
  if(action == "X"){
    break
  } else if (action == "L"){
    tic("finished LRT test with")
    cat(paste0("running deseq2 with LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
    
    dgelistl <- purrr::map2(summarized_exp, gsub("~", "", user_groups$full_model), function(x,y){dgel <- DGEList(counts=x$counts, group=x$sample_sheet[[y]]); return(dgel)})

    qlfl <- purrr::map2(summarized_exp, user_groups$full_model, run_edger)
    
    
    
    
    
    ddsl <- purrr::map2(user_groups$full_model, summarized_exp, function(x,y){dds <- DESeqDataSetFromMatrix(countData = y$counts,
                                                                                                            colData = y$sample_sheet,
                                                                                                            design = as.formula(x)); return(dds)})
    names(ddsl) <- user_groups$analysis_id
    
    # minimal filtering to remove outlier genes present in few cells ----------
    ddsl <- lapply(ddsl, function(x) x[cataract::meetsCondition(counts(x), gtet=10, nCells=10),])
    
    run_DESeq <- function(dds_obj, reduced_model){
      dds <- DESeq(dds_obj, reduced=as.formula(reduced_model), parallel = TRUE, BPPARAM = MulticoreParam(), betaPrior = FALSE, test="LRT")
      return(dds)
    }
    
    ddsl <- purrr::map2(ddsl, user_groups$reduced_model, run_DESeq)
    
    ddsRDS <- paste0(opt$out, user_groups$analysis_id, "_dds.rds")
    purrr::map2(ddsl, ddsRDS, saveRDS)
    ddsl <- lapply(ddsRDS, readRDS)
    
    get_res <- function(dds){
      res <- results(dds)
      res <- res[,-c(2,3)]
      baseMeanPerLvl <- sapply( levels(dds[[as.character(design(dds)[2])]]), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds[[as.character(design(dds)[2])]] == lvl] ) )
      res <- cbind(res, baseMeanPerLvl)
      genes <- data.frame(row.names = rownames(res), "gene_symbol" = lookup_genes(rownames(res)))
      res <- cbind(genes, res)
      return(res)
    }

    resl <- lapply(ddsl, get_res)
    
    print(paste("saving output csvs to", rescsv))
    rescsv <- paste0(opt$out, user_groups$analysis_id, ".csv")
    purrr::map2(resl, rescsv, write.csv)
    toc()
    
    
  } else if (action == "D"){
      
      check_file_exist <- function(destfile){
  
        if(file.exists(destfile)){
          dds <- readRDS(destfile)
        } else
          cat("preparing deseq2 input\n")
        }
  
      ddsRDS <- paste0(opt$out, user_groups$analysis_id, "_dds.rds")
  
      ddsl <- lapply(ddsRDS, check_file_exist)
  

      if(any(sapply(ddsl, is.null))){
        cat("running deseq2 for between group comparison")
        ddsl <- purrr::map2(user_groups$full_model, summarized_exp, function(x,y){dds <- DESeqDataSetFromMatrix(countData = y$counts,
                                                                                                                colData = y$sample_sheet,
                                                                                                                design = as.formula(x)); return(dds)})
        
        names(ddsl) <- user_groups$analysis_id
        
        # minimal filtering to remove outlier genes present in few cells ----------
        ddsl <- lapply(ddsl, function(x) x[cataract::meetsCondition(counts(x), gtet=10, nCells=10),])
        
        run_DESeq <- function(dds_obj){
          dds <- DESeq(dds_obj, parallel = TRUE, BPPARAM = MulticoreParam(), betaPrior = FALSE)
          return(dds)
        }
        
        ddsl <- purrr::map(ddsl, run_DESeq)
        ddsRDS <- paste0(opt$out, user_groups$analysis_id, "_dds.rds")
        purrr::map2(ddsl, ddsRDS, saveRDS)
        ddsl <- lapply(ddsRDS, readRDS)
        
      }
  
      cell_info <- cataract::safe_read(opt$cell_info, row.names = 1)
  
      results_models <- gsub("~", "", user_groups$full_model)
      results_models <- mapply(c, "full_model" = results_models, "analysis_id" = user_groups$analysis_id, SIMPLIFY = FALSE)
      
      diffex_comps <- purrr::map2(ddsl, results_models, cataract::make_all_deseq2_comparisons, cell_info) 

  } else if (action == "I"){
  browser()
  print(a)
  }
}

# plot normalized counts for a gene across deseq groups
# =========================================================
# gene_of_interest= "ENST00000371708" # enter the transcript id of the gene you are itnersted in!
# plotCounts(dds, gene=gene_of_interest, intgroup=c("treatment_group", "day"))
# 
# 
# report <- DESeq2Report(ddsMF, 'DESeq2Report-sunhye_all_treat_groups', 'treatment_group',
#                        outdir = 'DESeq2-Report')
#
# descriptions <- lapply(test, function(x)x@elementMetadata@listData[["description"]][[2]] )
#
# comparison_names <- lapply(strsplit(unlist(descriptions), " "),
#                            function(x) paste(x[c(6,7,8)], collapse="_"))
#
# test <- report_deseq2_results(ddsMF, "day")
#
# test0 <- filter_deseq_results(test)
#
# day_deseq_results <- purrr::map(test0, clean_deseq_trs_ediff)
#
# names(day_deseq_results) <- comparison_names
#
# purrr::map2(day_deseq_results, paste0(names(day_deseq_results), "_deseq2_results.csv"), write.table)
#
# saveRDS(day_deseq_results, "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/day_deseq_results.rds")
#
#
# library(IHW)
# resIHW <- results(dds, filterFun=ihw)
# summary(resIHW)
#
# sum(resIHW$padj < 0.1, na.rm=TRUE)
# metadata(resIHW)$ihwResult
#
# resultsNames(dds)
#
# d <- plotCounts(dds, gene=which.min(results_table$padj), intgroup="treatment_group",
#                 returnData=TRUE)



