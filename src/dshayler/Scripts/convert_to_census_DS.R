#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<=1) {
  stop("At least two arguments must be supplied 1) input file 2) sample sheet file", call.=FALSE)
} 

# load required libraries
library(reshape2)
library(purrr)
library(dplyr)
library(tibble)
require(gtools)
library(monocle)

EXPRESSION_PATH = "single_cell_pipeline/output/FACS_20171031_dshayler_H_sapiens_output/stringtie_transcripts.tpm.csv" #this points to Seq_2 gene expression file
EXPRESSION_PATH = c(EXPRESSION_PATH, "single_cell_pipeline/output/FACS_20170407_dshayler_H_sapiens_output/stringtie_transcripts.tpm.csv") #this points to Seq_1 gene expression file

jt_complete = map(EXPRESSION_PATH, read.table, header = TRUE, sep="\t") %>% 
  reduce(merge, by="TRANSCRIPT_ID") %>% 
  remove_rownames %>% 
  column_to_rownames(var="TRANSCRIPT_ID")

census_output_file <- gsub(".csv", "_census_matrix.csv", EXPRESSION_PATH[[1]])

GROUP_PATH = "single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20171031_dshayler_H_sapiens/Sample_number_attributes_103117.csv"
GROUP_PATH = c(GROUP_PATH, "single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20170407_dshayler_H_sapiens/Dominic_Cell_Age_and_Sort_Method_042017.csv")

cell_groups = map_dfr(GROUP_PATH, read.table, header = TRUE, sep = ",")
rownames(cell_groups) <- cell_groups[,1]

# read in data and set up objects

stringtie_transcripts_matrix <- jt_complete
colClean <- function(x){ colnames(x) <- gsub("_.*", "", colnames(x)); x }
stringtie_transcripts_matrix <- colClean(stringtie_transcripts_matrix)
stringtie_transcripts_matrix <- stringtie_transcripts_matrix[ , mixedsort(colnames(stringtie_transcripts_matrix))]
sample_sheet <- cell_groups
transcripts_annotation <- data.frame("transcripts" = rownames(stringtie_transcripts_matrix), row.names = rownames(stringtie_transcripts_matrix))



# check if rownames of sample_sheet and colnames of stringtie_transcript_matrix match
if(any(rownames(sample_sheet) != colnames(stringtie_transcripts_matrix))){
  stop("Check that format of sample sheet is correct and that all samples were included in analysis")
}

pd <- new("AnnotatedDataFrame", data=sample_sheet)
fd <- new("AnnotatedDataFrame", data=transcripts_annotation)


# First create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(as.matrix(stringtie_transcripts_matrix),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=0.1,
                       expressionFamily=tobit(Lower=0.1))
# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM)
# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# # keep genes that are expressed in at least 5 cells
# min_expression <- 0.1
# HSMM <- detectGenes(HSMM, min_expr=min_expression)
# expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 5))

#print output of census to csv prior to monocle workflow
write.table(as.matrix(exprs(HSMM)), census_output_file, sep="\t")

