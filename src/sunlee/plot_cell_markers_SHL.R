#!/usr/bin/Rscript


# load required libraries -------------------------------------------------

library(biomaRt)
library(shiny)
library(heatmaply)
library(shinyHeatmaply)

ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
t_to_g_association <- readRDS("~/Homo_sapiens/grch38_tran/t_to_g_assocation.rds")

cone_markers <- c("OPN1LW", "RXRG", "CRX", "ARR3", "RORB", "GTF2IRD1", "GNAT2", "THRB", "MDM2", "MYCN")
rod_markers <- c("NRL", "NXNL1")
ama_hoz_markers <- c("PROX1", "PAX6", "NES")
prog_mul_markers <- c("NES", "RLBP1", "SOX2")

# load required functions -------------------------------------------------

gene_sym_from_trsid <- function(matrix, mart){
  
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = rownames(matrix), mart = mart)
  biomaRt_result <- as.data.frame(biomaRt_result, stringsAsFactors = FALSE)
  gene = biomaRt_result$hgnc_symbol
  ens_trans = biomaRt_result$ensembl_transcript_id
  names(ens_trans) = gene
  assign("t_to_g_association",ens_trans,envir = .GlobalEnv)
}

#finds genename in ensemble<->gene association tables
getGene <- function(gene, association = t_to_g_association) {
  gene = paste0("^",gene,"$")
  grep_results <- grep(gene, names(association))
  transcript_ids <- unname(association[grep_results])
  #results <- (as.character(association)[lapply(names(t_to_g_association), function(x) which(gene%in%x))])
  return(transcript_ids)
}

SHL_0407 <- read.table("~/single_cell_pipeline/results/sunhye/sunhye_census_matrix.csv", sep = "\t", header = TRUE)

SHL_1031 <- read.table("~/single_cell_pipeline/output/FACS_20171031_sunlee_H_sapiens_output/transcripts.tpm_census_matrix.csv", sep = "\t", header = TRUE)

# cone_df <- purrr::map_df(cone_markers, function(x){setNames(colSums(SHL_1031[getGene(x),]), x)})

return_counts <- function(gene_vec, expr_mat){

  expr_list <- lapply(gene_vec, function(x){expr_mat[getGene(x),]})
  names(expr_list) <- gene_vec
  gene_count_df <- purrr::map_dfr(expr_list, colSums)
  gene_count_df <- as.data.frame(gene_count_df)
  row.names(gene_count_df) <- colnames(expr_mat)
  return(gene_count_df)
}

marker_list <- list("cone_markers" = cone_markers, "rod_markers" = rod_markers, "ama_hoz_markers" = ama_hoz_markers, "prog_mul_markers" = prog_mul_markers)

gene_count_df <- purrr::map(marker_list, return_counts, SHL_0407)

gene_count_mat <- lapply(gene_count_df, data.matrix)

save(SHL_0407_gene_count_mat, file = "~/single_cell_pipeline/results/sunhye/gene_count_mat_0407.rda")
save(SHL_1031_gene_count_mat, file = "~/single_cell_pipeline/results/sunhye/gene_count_mat_1031.rda")


# create separate csv files for each dataframe in a list ------------------
sapply(names(SHL_1031_gene_count_mat), 
       function (x) write.table(SHL_1031_gene_count_mat[[x]], file=paste(x, "txt", sep=".") )   )


