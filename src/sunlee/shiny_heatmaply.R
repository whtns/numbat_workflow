#!/usr/bin/Rscript


# load required libraries -------------------------------------------------

library(shiny)
library(heatmaply)
library(shinyHeatmaply)
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

lookup_genes <- function(txname){
  
  txs <- transcripts(edb, filter = TxIdFilter(txname),
                     columns = c("symbol"))
  return(txs$symbol)
  
}

lookup_transcripts <- function(genename){
  
  txs <- transcripts(edb, filter = GenenameFilter(genename),
                     columns = c("tx_id"))
  return(txs$tx_id)
  
}



# load census matrix ------------------------------------------------------

census_matrix <- "~/single_cell_pipeline/output/dshayler/20170407_20171031_combined_census_matrix.tpm.csv"
census_matrix <- read.table(census_matrix, sep = "\t", header = TRUE)

runApp(system.file("shinyapp", package = "shinyHeatmaply"))

mat_path <- "~/single_cell_pipeline/results/sunhye/gene_count_mat_0407.rda"

reorder_mat <- function(mat_path, marker_class, sort_gene){
  browser()
  x <- load(mat_path)
  sort_mat <- get(x)
  rm(x)
  
  sort_mat[[marker_class]] <- sort_mat[[marker_class]][order(sort_mat[[marker_class]][,sort_gene], decreasing = TRUE),] 
  save(sort_mat, file = mat_path)
  return(sort_mat)
}


ddsRDS <- "~/single_cell_pipeline/output/merged_analyses/DESEQ2/dds_cluster.rds"
saveRDS(dds, ddsRDS)
dds <- readRDS(ddsRDS)

contrasts <- combn(levels(sample_sheet$cluster), 2)
contrasts <- purrr::map2(contrasts[1,], contrasts[2,], c)

# resl <- purrr::map(contrasts[1:2], function(x) results(dds, contrast = c("cluster", x[1], x[2])))
# 
# resl <- lapply(resl, as.data.frame)
# resl <- lapply(resl, function(x)na.omit(x[(x$padj < 0.05),]))
# 
# res_genes <- lookup_genes(rownames(resl[[1]]))
# 
# resl <- lapply(resl, function(x)cbind("gene_id" = res_genes, x))


res <- results(dds)
res <- as.data.frame(res)
res <- na.omit(res[(res$padj < 0.05),])


filt_census_matrix <- census_matrix[rownames(census_matrix) %in% rownames(res),]

filt_census_matrix <- t(apply(filt_census_matrix, 1, heatmaply::normalize))

filt_census_matrix <- filt_census_matrix[,mixedsort(labels(col_dend))]



jt.l_tmp <- readRDS("~/single_cell_pipeline/src/dominic/Scripts/jt.l_tmp")

col_dend  <- grx %>%
  color_branches(k=9) %>% 
  ladderize

heatmaply(assay(vsd)[ topVarGenes, ], scale = "row", seriate = "none", dendrogram = "column", Colv = col_dend, hclust_method = "complete", file = "./heatmaply_plot.html")

# use the log transform on the data set
vsd <- vst(dds, blind=F)

#cell by cell heatmap
sampleDists <- dist( t( assay(vsd) ) )
sampleDistMatrix <- as.matrix( sampleDists )

colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)


# gene by cell heatmap ----------------------------------------------------


topVarGenes <- head( order( rowVars( assay(vsd) ), decreasing=TRUE ), 100 )
heatmap.2( assay(vsd)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


ramp <- 1:3/3
cols <- c(rgb(ramp, 0, 0),
          rgb(0, ramp, 0),
          rgb(0, 0, ramp),
          rgb(ramp, 0, ramp))
print( plotPCA( vsd, intgroup = c( "wet", "d"), rycol=cols ) )

heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))


library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1),
           ylim = c(0,2.5))
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))

# scratchpad --------------------------------------------------------------



topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=T),100)
matrix <- assay(vsd)[ topVarGenes ]
matrix <- matrix - rowMeans(matrix)

# select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("ConditionA","ConditionB")])
pheatmap(matrix, annotation_col=annotation_data)

topVarGenes <- head(order(-rowVars(filt_census_matrix)),35)


browseURL("./heatmaply_plot.html")
