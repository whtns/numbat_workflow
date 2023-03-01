#!/usr/bin/Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(biomaRt)
library(tibble)
library(gtools)

gene_sym_from_trsid <- function(x, grp_comp=grp_comp){
  browser()
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = x$`ensembl_transcript_id`, mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  output <- setDT(as.data.frame(x), keep.rownames = TRUE)[]
  output <- rename(output, ensembl_transcript_id = rn)
  output <-  merge(as.data.frame(biomaRt_result), output, by = "ensembl_transcript_id") 
  return(output)
}

tophat_alignment_file <- "./stats/cufflinks_stats_2015_2016/alignment_statistics_new.csv"

old_cell_group_file <- "./stats/cufflinks_stats_2015_2016/cell_technology_new_22-81.csv"
old_alignment_stats_file = "./stats/SHL_201607_rerun_hisat2_alignment_stats.csv"
old_stringtie_gene_counts_file = "./stats/cufflinks_stats_2015_2016/gene_count_matrix.csv"



old_cell_groups <- read.table(old_cell_group_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
old_reads_table <- read.table(old_alignment_stats_file, header = TRUE, sep = "\t")
old_stringtie_gene_counts <- read.table(old_stringtie_gene_counts_file, header = TRUE, sep=",", row.names = 1)
old_stringtie_gene_counts <- stringtie_gene_counts[,mixedsort(colnames(stringtie_gene_counts))]

new_cell_group_file <- "./stats/SAMPLE_SHEET_hiseq_SHL_04172017.csv"
new_alignment_stats_file = "./stats/cufflinks_stats_2015_2016/hiseq_alignment_stats.csv"
new_stringtie_gene_counts_file = "./output/gene_count_matrix.csv"

new_cell_groups <- read.table(new_cell_group_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
new_reads_table <- read.table(new_alignment_stats_file, header = TRUE, sep = "\t")
new_stringtie_gene_counts <- read.table(new_stringtie_gene_counts_file, header = TRUE, sep=",", row.names = 1)
new_stringtie_gene_counts <- stringtie_gene_counts[,mixedsort(colnames(stringtie_gene_counts))]


new_cell_groups <- new_cell_groups %>%
  mutate(cell_name = gsub("^X", "Hu_", Sample_ID)) %>%
  dplyr::select(cell_name) %>%
  group_by(cell_name) %>%
  mutate(chemistry = "HiSeq_w_Smarter_v4")

tophat_alignment_stats <- read.table(tophat_alignment_file, header = TRUE) 
tophat_gene_counts_filenames <- list.files(path="./sc_cone_devel/sc_cone_devel_H_sapiens/C1_20160718_SHL_H_sapiens/output_cufflinks_MT/", 
                                           pattern="genes.fpkm_tracking", full.names  = TRUE, recursive = TRUE)
tophat_gene_counts_list <- lapply(tophat_gene_counts_filenames, function(x)read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill=TRUE))
tophat_gene_counts_names <-substr(tophat_gene_counts_filenames,91,92) 

names(tophat_gene_counts_list) <- tophat_gene_counts_names
tidy_tophat <- tophat_gene_counts_list%>%
  rbindlist(use.names=TRUE, fill=TRUE, idcol="cell_name") %>%
  mutate(cell_name = paste0("Hu_", cell_name)) %>%
  inner_join(old_cell_groups, by = "cell_name") %>%
  mutate(chemistry = paste0("tophat_", chemistry))

# filter and merge alignent tables ----------------------------------------

old_reads_table <- old_reads_table %>%
  mutate(cell_name = paste0("Hu_", X)) %>%
  inner_join(old_cell_groups, by = "cell_name") %>%
  dplyr::select(cell_name, chemistry, PAIR.AL.C.1) %>%
  mutate(PAIR.AL.C.1 = 2 * PAIR.AL.C.1)

new_reads_table <- new_reads_table %>%
  mutate(cell_name = paste0("Hu_", X)) %>%
  inner_join(new_cell_groups, by = "cell_name") %>%
  dplyr::select(cell_name, chemistry, PAIR.AL.C.1) %>%
  mutate(PAIR.AL.C.1 = 2 * PAIR.AL.C.1)

tophat_alignment_stats <- tophat_alignment_stats %>%
  mutate(chemistry = "tophat_alignment")
  dplyr::select(cell_name, chemistry, concordant)

chemistries <- c("HiSeq_w_Smarter_v4", "SMARTer_V4", "tophat_alignment")
summ_reads_table <- bind_rows(old_reads_table, new_reads_table) %>%
  rename(cell = cell_name, concordant = PAIR.AL.C.1) %>%
  bind_rows(tophat_alignment_stats) %>%
  filter(chemistry %in% chemistries) %>%
  rename(concordant_mapped_reads = concordant)


# filter and merge gene count tables --------------------------------------
tophat_gene_counts <- tidy_tophat %>%
  filter(FPKM > 0) %>%
  group_by(cell_name, chemistry) %>%
  summarise(genes_detected = n()) 

old_gene_counts <- old_stringtie_gene_counts %>% 
  rownames_to_column(var = "ensembl_transcript_id") %>%
  gather(cell_name, raw_count, -ensembl_transcript_id) %>%
  inner_join(old_cell_groups, by = "cell_name") %>%
  group_by(cell_name, chemistry) %>%
  filter(raw_count > 0) %>%
  summarise(genes_detected = n())

new_gene_counts <- new_stringtie_gene_counts %>% 
  rownames_to_column(var = "ensembl_transcript_id") %>%
  gather(cell_name, raw_count, -ensembl_transcript_id) %>%
  mutate(cell_name = gsub("^X", "Hu_", cell_name)) %>%
  inner_join(new_cell_groups, by = "cell_name") %>%
  group_by(cell_name, chemistry) %>%
  filter(raw_count > 0) %>%
  summarise(genes_detected = n())

chemistries <- c("HiSeq_w_Smarter_v4", "SMARTer_V4", "tophat_SMARTer_V4")
summ_gene_counts <- bind_rows(old_gene_counts, new_gene_counts) %>%
  bind_rows(tophat_gene_counts) %>%
  filter(chemistry %in% chemistries) 

# do plotting -------------------------------------------------------------

plot_list = list()
pdf("./stats/tophat_stats_2015_2016/box_plots_method_compare.pdf", onefile = TRUE)
par(mfrow = c(2,1))

ggplot(summ_reads_table, aes(x=chemistry, y = concordant_mapped_reads)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggplot(summ_gene_counts, aes(x=chemistry, y=genes_detected)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dev.off()







