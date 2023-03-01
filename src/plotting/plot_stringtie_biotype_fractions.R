#!/usr/bin/Rscript


# load required libraries -------------------------------------------------

library(dplyr)
library(cowplot)


# load datasets -----------------------------------------------------------

stringtie_trs_matrix <- read.table("../output/transcript_count_matrix.csv", sep =",", header = T, row.names = 1)
stringtie_gene_matrix <-  read.table("../output/gene_count_matrix.csv", sep = ",", header = T, row.names = 1)

transcript_biotypes <- read.table("../../Homo_sapiens/grch38_tran/transcript_biotypes.csv", stringsAsFactors = FALSE, header = TRUE)
gene_biotypes <- read.table("../../Homo_sapiens/grch38_tran/gene_biotypes.csv", stringsAsFactors = FALSE, header = TRUE)


# load required functions -------------------------------------------------

meetsCondition = function(mat, gtet = 0, nCells=0){
  condition_mat = ifelse(mat>=gtet,1,0)
  meets_condition = apply(condition_mat,1,sum) >= nCells
  return(meets_condition)
}

get_biotype_fractions <- function(expression_matrix, biotype_table, id_type){
  
  raw_biotype_fractions <- expression_matrix %>% 
    rownames_to_column(var = id_type) %>% 
    mutate(varsum = rowSums(dplyr::select(., contains("X")))) %>% 
    left_join(biotype_table, by = id_type) %>% 
    mutate(coding_status = ifelse(grepl("protein_coding", gene_biotype), "coding", "noncoding"))
  
  meets_condition = meetsCondition(raw_biotype_fractions,gtet=0.1,nCells = 10)
  raw_biotype_fractions = as.data.frame(raw_biotype_fractions[meets_condition,])
  
  give.n <- function(x){
    return(c(y = min(x)/2, label = length(x)))
  }
  
  coding_binary_plot <- ggplot(raw_biotype_fractions, aes(varsum, color=coding_status)) + geom_freqpoly()
  print(coding_binary_plot)
  
  raw_biotype_plot <- ggplot(raw_biotype_fractions, aes(gene_biotype, varsum)) + geom_boxplot() +
    coord_flip() + stat_summary(fun.data = give.n, geom = "text", cex = 2.0, color = "red") +
    theme(axis.text = element_text(size = 5)) + labs(title = id_type)
  print(raw_biotype_plot)
  
  filt_biotype_fractions <- dplyr::filter(raw_biotype_fractions, varsum < 30000)
  filt_biotype_plot <- raw_biotype_plot %+% filt_biotype_fractions
  print(filt_biotype_plot)
  
  combined_plot <- plot_grid(coding_binary_plot, raw_biotype_plot, filt_biotype_plot, align = "hv", nrow =3) + 
    guides(colour = guide_legend(nrow = 2)) 
  print(combined_plot)
  ggsave(paste0("../results/sunhye/rna_biotype_distribution_shl_", id_type, ".pdf"), combined_plot, width = 6, height = 8)
  
}


# run plotting ------------------------------------------------------------

get_biotype_fractions(stringtie_trs_matrix, transcript_biotypes, "ensembl_transcript_id")

get_biotype_fractions(stringtie_gene_matrix, gene_biotypes, "ensembl_gene_id")
