#!/usr/bin/Rscript


# Caution!  ---------------------------------------------------------------

# known bug in scde error fitting, see https://groups.google.com/forum/#!topic/singlecellstats/rbFUTOQ9wu4
#   need to install prior version of package 'flexmix'


# load required libraries ----------------------------------------------
library(scde)
library(tidyverse)
library(gtools)
library(biomaRt)
library(data.table)
library(regionReport)
library(dplyr)
library(cowplot)

# load required functions -------------------------------------------------

query_biomart_trs <- function(x, query_col, ...){
 
  query_col <- deparse(substitute(query_col))
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = x[[query_col]], mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  output <-  left_join(as.data.frame(biomaRt_result), x, by = c("ensembl_transcript_id" = query_col)) 
  }

query_biomart_genes <- function(x, query_col, ...){
  
  query_col <- deparse(substitute(query_col))
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_gene_id', values = x[[query_col]], mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  output <-  left_join(as.data.frame(biomaRt_result), x, by = c("ensembl_gene_id" = query_col)) 
}

stringtie_trs_matrix <- read.table("../output/transcript_count_matrix.csv", sep =",", header = T, row.names = 1)
stringtie_gene_matrix <-  read.table("../output/gene_count_matrix.csv", sep = ",", header = T, row.names = 1)

branch_GROUP_PATH = "../data/Sunhye_cell_division_day_treatment_and_branch.csv"

#Make the gene-wise matrix
coldata <- read.table(branch_GROUP_PATH, sep=",", header = TRUE, stringsAsFactors = FALSE)

prep_scde_input <- function(expr_matrix, cell_data, comparison, feature){
  browser()
  expr_matrix <- expr_matrix[,colnames(expr_matrix)%in%cell_data$cell_id] 
  
  expr_matrix <- expr_matrix[ , mixedsort(colnames(expr_matrix))]
  
  if(all(rownames(cell_data) == colnames(expr_matrix))){
    scde_input <- list("expr_matrix" = expr_matrix, "cell_data" = cell_data)
  } else{
    stop("cells in expression matrix and cell data file must match; recheck input")
  }
  

  # begin scde analysis -----------------------------------------------------
  
  group_factor <- factor(cell_data[[comparison]])
  names(group_factor) <- colnames(expr_matrix)  
  table(group_factor)
    
  cd <- clean.counts(expr_matrix, min.lib.size = 1000, min.reads = 1, min.detected = 1)
  
  # calculate models
  o.ifm <- scde.error.models(counts = cd, groups = group_factor, n.cores = 7, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
  
  # filter out cells that don't show positive correlation with
  # the expected expression magnitudes (very poor fits)
  valid.cells <- o.ifm$corr.a > 0
  table(valid.cells)
  o.ifm <- o.ifm[valid.cells, ]
  
  # estimate gene expression prior
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
  
  scde_input <- list("group_factor" = group_factor, "cd" = cd, "o.ifm" = o.ifm, "o.prior" = o.prior)
  saveRDS(scde_input, paste0("../results/sunhye/diff_expression/diffex_by", feature, "/", comparison, "_", feature, "_stringtie_sunhye.rds"))
  return(scde_input)
}

branch_trs_scde_input <- prep_scde_input(stringtie_trs_matrix, coldata, "branch", "trs")

branch_gene_scde_input <- prep_scde_input(stringtie_gene_matrix, coldata, "branch", "gene")

sub_branch_trs_scde_input <- prep_scde_input(stringtie_trs_matrix, coldata, "sub_branch", "trs")

sub_branch_gene_scde_input <- prep_scde_input(stringtie_gene_matrix, coldata, "sub_branch", "gene")

parallel_gene_scde <- function(scde_input, comparison){
  browser()
  run_gene_scde <- function(comparison_groups, scde_input){
    # define two groups of cells
    browser()
    groups <- scde_input$group_factor
    names(groups) <- row.names(scde_input$o.ifm)
    groups[!groups %in% comparison_groups] <- NA
    groups <- factor(groups)
    ediff <- scde.expression.difference(scde_input$o.ifm, scde_input$cd, scde_input$o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  7, verbose  =  1)
    group_levels <-  paste(levels(groups), collapse = "_")
    saveRDS(ediff, paste0("../results/sunhye/diff_expression/diffex_by_gene/gene_branch_", group_levels, "_diffex.rds"))
    
    clean_gene_ediff <-function(ediff){
      group_ediff <- ediff %>% 
        rownames_to_column(var = "ensembl_gene_id") %>% 
        mutate(p_val=2*pnorm(-abs(Z))) %>% 
        arrange(p_val) %>% 
        dplyr::select(ensembl_gene_id, mle, p_val)
      
      group_ediff <- query_biomart_genes(group_ediff, ensembl_gene_id)
    }
    
    ediff <- clean_gene_ediff(ediff)
    
    ediff <- arrange(ediff) %>% 
      arrange(desc(mle))
    
  }
  
  if (comparison == "branch"){
    diffex_comparisons <- combn(levels(scde_input$group_factor), 2)
    diffex_comparisons <- map2(diffex_comparisons[1,], diffex_comparisons[2,], c)
    diffex_names <- map(diffex_comparisons, paste, collapse="_vs_")
    diffex_paths <- map(diffex_names, function(x){paste0("../results/sunhye/diff_expression/diffex_by_gene/branch_", x, "_stringtie_gene_expression_values.csv")})
  } else if (comparison == "sub_branch"){
    
    groups <- scde_input$group_factor
    groups[groups %in% c("a1", "a2")] <- NA
    groups <- factor(groups)
    
    diffex_comparisons <- map(levels(groups), function(x){c("a1", x)})
    diffex_names <- map(diffex_comparisons, paste, collapse="_vs_")
    diffex_paths <- map(diffex_names, function(x){paste0("../results/sunhye/diff_expression/diffex_by_gene/sub_branch_", x, "_stringtie_gene_expression_values.csv")})
  } else {
    stop("invalid input")
  }
  
  gene_comparison_ediff <- map(diffex_comparisons, run_gene_scde, scde_input)
  gene_comparison_ediff <- setNames(gene_comparison_ediff, diffex_names)
  
  map2(gene_comparison_ediff, diffex_paths, write.table)
  
}

parallel_trs_scde <- function(scde_input, comparison){
  browser()
  run_gene_scde <- function(comparison_groups, scde_input){
    browser()
    # define two groups of cells
    groups <- scde_input$group_factor
    names(groups) <- row.names(scde_input$o.ifm)
    groups[!groups %in% comparison_groups] <- NA
    groups <- factor(groups)
    ediff <- scde.expression.difference(scde_input$o.ifm, scde_input$cd, scde_input$o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  7, verbose  =  1)
    group_levels <-  paste(levels(groups), collapse = "_")
    saveRDS(ediff, paste0("../results/sunhye/diff_expression/diffex_by_trs/trs_branch_", group_levels, "_diffex.rds"))
    
    clean_trs_ediff <-function(ediff){
      group_ediff <- ediff %>% 
        rownames_to_column(var = "ensembl_transcript_id") %>% 
        mutate(p_val=2*pnorm(-abs(Z))) %>% 
        arrange(p_val) %>% 
        dplyr::select(ensembl_transcript_id, mle, p_val)
      
      group_ediff <- query_biomart_transcripts(group_ediff, ensembl_transcript_id)
    }
    
    ediff <- clean_trs_ediff(ediff)
    
    ediff <- arrange(ediff) %>% 
      arrange(desc(mle))
    
  }
  
  if (comparison == "branch"){
    diffex_comparisons <- combn(levels(scde_input$group_factor), 2)
    diffex_comparisons <- map2(out[1,], out[2,], c)
    diffex_names <- map(diffex_comparisons, paste, collapse="_vs_")
    diffex_paths <- map(diffex_names, function(x){paste0("../results/sunhye/diff_expression/diffex_by_trs/branch_", x, "_stringtie_trs_expression_values.csv")})
  } else if (comparison == "sub_branch"){
    
    groups <- scde_input$group_factor
    groups[groups %in% c("a1", "a2")] <- NA
    groups <- factor(groups)
    
    diffex_comparisons <- map(levels(groups), function(x){c("a1", x)})
    diffex_names <- map(diffex_comparisons, paste, collapse="_vs_")
    diffex_paths <- map(diffex_names, function(x){paste0("../results/sunhye/diff_expression/diffex_by_trs/sub_branch_", x, "_stringtie_trs_expression_values.csv")})
  } else{
    stop("invalid input")
  }
  
  trs_comparison_ediff <- map(diffex_comparisons, run_trs_scde, scde_input)
  trs_comparison_ediff <- setNames(trs_comparison_ediff, diffex_names)
  
  map2(trs_comparison_ediff, diffex_paths, write.table)
  
}


# run gene scde with comparison between branches --------------------------

scde_input <- readRDS("../results/sunhye/diff_expression/branch_gene_stringtie_sunhye.rds")
parallel_gene_scde(scde_input, "branch")


# run gene scde with comparison between sub_branches ----------------------

scde_input <- readRDS("../results/sunhye/diff_expression/sub_branch_gene_stringtie_sunhye.rds")
parallel_gene_scde(scde_input, "sub_branch")

# run trs scde with comparison between branches --------------------------

scde_input <- readRDS("../results/sunhye/diff_expression/branch_gene_stringtie_sunhye.rds")
parallel_trs_scde(scde_input, "branch")


# run trs scde with comparison between sub_branches ----------------------

scde_input <- readRDS("../results/sunhye/diff_expression/sub_branch_gene_stringtie_sunhye.rds")
parallel_trs_scde(scde_input, "sub_branch")



# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])


scde.test.gene.expression.difference("PROZ", models = o.ifm, counts = cd, prior = o.prior)


# tidyverse comparisions --------------------------------------------------

testa_b <- rownames_to_column(ediff_a_b, "gene_id") %>%
  top_n(1000, mle) %>% 
  arrange(desc(mle))

testa_c <- rownames_to_column(ediff_a_c, "gene_id") %>% 
  top_n(1000, mle) %>%
  rename(mle_ac = mle) %>% 
  select(gene_id, mle_ac) 

test <- full_join(testa_b, testa_c) 
  gather("mle", "diff_comp", 2:3)

ggplot(test, aes(x = gene_id))
setdiff(testa_b$gene_id, testa_c$gene_id)
intersect(testa_b$gene_id, testa_c$gene_id)


