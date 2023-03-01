#!/usr/bin/Rscript

#load required libraries and functions
#=====================================
library(Rsubread)
library(DESeq2)
library(tidyverse)
library(gtools)
library(biomaRt)
library(data.table)
library(regionReport)

source("/home/skevin/TOOLS/maxwellbay-waterfall-9d035d9cf75a/src/HyperGeoTerms.R")

gene_sym_from_trsid <- function(x, grp_comp=grp_comp){
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = rownames(x), mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  output <- setDT(as.data.frame(x), keep.rownames = TRUE)[]
  output <- rename(output, ensembl_transcript_id = rn)
  output <-  merge(as.data.frame(biomaRt_result), output, by = "ensembl_transcript_id") 
  filename = paste0(grp_comp,"_wilcox_raw_count",".csv")
  write.table(output,filename, sep="\t", quote=FALSE, row.names = FALSE)
}
#Generate RNA-seq matrix
#Set parameters
#=================================================

GTF="/media/thor/storage/os_exosome_pipeline/gencode.v26lift37.annotation.gtf"
#stringtie_transcripts_file <- "/dataVolume/storage/single_cell_pipeline/stats/Sunhye_stringtie.tpm.csv"
sample_sheet_file <- "/dataVolume/storage/single_cell_pipeline/stats/SAMPLE_SHEET_hiseq_SHL_04172017.csv"
transcripts_attr_file <- "/dataVolume/storage/single_cell_pipeline/stats/Sunhye_stringtie_transcript_annotation.csv"

# read in data and set up objects
#stringtie_transcripts_matrix <- read.table(stringtie_transcripts_file, header=T, sep=",", row.names=1)
#stringtie_transcripts_matrix <- stringtie_transcripts_matrix[ , mixedsort(colnames(stringtie_transcripts_matrix))]
sample_sheet <- read.table(sample_sheet_file, header=T, sep=",", row.names=1)
transcript_ann <- read.table(transcripts_attr_file, row.names=1, header=T, sep=",")

#Make the gene-wise matrix
coldata <- sample_sheet
counts <- read.table("/dataVolume/storage/single_cell_pipeline/output/transcript_count_matrix.csv", sep =",", header = TRUE, row.names = 1)




# analysis with Waterfall::signatureAll -----------------------------------

counts_without_controls = counts[,coldata[,"treatment_group"] != "shCtrl"]
coldata_reduced <- coldata[colnames(counts_without_controls),]
signature_vector <- rownames(coldata_reduced)
names(signature_vector) <- coldata_reduced[,2]
signature_vector <- signature_vector[names(signature_vector) == "sh733"]
wilcox_sig_matrix <- signatureAll(counts_without_controls, signature_vector)
wilcox_adjust <- p.adjust(wilcox_sig_matrix$pval, "BH")
wilcox_2 <- wilcox_sig_matrix %>%
  add_column(p.adjust = wilcox_adjust, .after = 2) %>%
  filter(meanA<Inf, log2FC<Inf, !is.na(log2FC)) 
wilcox_2 <- as.data.frame(wilcox_2)
rownames(wilcox_2) <- wilcox_2$name
gene_sym_from_trsid(wilcox_2, "733_over_737")
  




# continue with deseq2 analysis -------------------------------------------
counts <- counts[ , mixedsort(colnames(counts))]
colnames(counts) <- rownames(coldata)
all(rownames(coldata) == colnames(counts))


dds_parallel <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ treatment_group)


dds_parallel <- dds_parallel[ rowSums(counts(dds_parallel)) > 1, ]

dds_parallel$`treatment_group` <- factor(dds_parallel$`treatment_group`, levels=c("shCtrl","sh733", "sh737"))

dds_parallel <- DESeq(dds_parallel, parallel = TRUE, betaPrior = FALSE)


# if no changes made to code above, enter command:
# dds <- readRDS("/dataVolume/storage/single_cell_pipeline/sunhye_deseq_dds.rds")
# ===============================

# plot normalized counts for a gene across deseq groups
# =========================================================
gene_of_interest= "ENST00000371708" # enter the transcript id of the gene you are itnersted in!
plotCounts(dds, gene=gene_of_interest, intgroup=c("treatment_group", "day"))


report <- DESeq2Report(dds, 'DESeq2Report-sunhye_all_treat_groups', 'treatment_group',
                       outdir = 'DESeq2-Report')


# to generate a table of results from differential expression analysis
#======================================================================
results_parallel <- results(dds_parallel, contrast = c("treatment_group", "sh733", "sh737"), parallel = TRUE)

#multi-factor design
#=============================
ddsMF <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ treatment_group + day)

ddsMF <- DESeq(ddsMF, parallel = TRUE)


resOrdered <- results_table[order(results_table$padj),]
resOrdered_l2fc <- results_table[order(results_table$log2FoldChange),]
summary(results_table)
sum(results_table$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE)

library(IHW)
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)

sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

resultsNames(dds)

d <- plotCounts(dds, gene=which.min(results_table$padj), intgroup="treatment_group", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=treatment_group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

plot(metadata(results_table)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(results_table)$lo.fit, col="red")
abline(v=metadata(results_table)$filterTheta)


