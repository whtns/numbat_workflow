#!/usr/bin/Rscript --vanilla --slave

# script provided by Matthew Thornton of Grubbs lab at USC HSC

today <- Sys.Date()
dt <- format(today, format="%Y%m%d")

library(Biobase)
library(org.Hs.eg.db)
library(GO.db)
library(annotate)
library(GOstats)
library(genefilter)
library(xtable)
library(Rgraphviz)
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86
library(topGO)

lookup_entrez_from_ensembl_geneid <- function(geneid){
  entrez <- ensembldb::genes(edb, filter = GeneIdFilter(geneid), columns = c("entrezid"))
  # some ensembl geneids have multiple entrez ids; drop all but first
  return(entrez)
}

## Change to data directory. This script is for Mouse, change Mm to Hs for Human.

diffex_file <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/EDGER/1_edger_lrt.csv"

## input reads for filtering for GOstats NCBI IDs are best for rownames, nicknamed "EntezID". Can't be duplicated
diffex <- read.csv(diffex_file, row.names = 1)

q_data <- diffex[diffex$FDR < 0.05, ]
# up_q <- q_data[q_data$logFC > 0, ]
# dn_q <- q_data[q_data$logFC < 0, ]

prep_topgo_data <- function(diffex_features){
  browser()
  
  entrezid <- lookup_entrez_from_ensembl_geneid(rownames(diffex_features))
  
  geneID2GO <- entrezid$entrezid
  
  GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  
  diffex$entrezid <- entrezid
  
  diffex <- diffex[!is.na(diffex$entrezid),]
  rownames(diffex) <- diffex$entrezid
  
  ## Convert to ExpressionSet
  diffex_features <- new("ExpressionSet", exprs=as.matrix(diffex_features))
  
  # filter by posession of a valid GO term
  haveGO <- sapply(mget(featureNames(diffex_features), org.Hs.egGO), function(x) { if (length(x) == 1 && is.na(x)) FALSE else TRUE})
  diffex_features <- diffex_features[haveGO, ]
  
  ## Format the results for using with GOstats
  diffex_features <- unique(unlist(featureNames(diffex_features)))
  
  return(diffex_features)
  
}

prep_gostats_data <- function(diffex_features){
  browser()
  ## Convert to ExpressionSet
  entrezid <- lookup_entrez_from_ensembl_geneid(rownames(diffex_features))
  entrezid <- lapply(entrezid$entrezid, `[[`, 1)
  
  diffex_features$entrezid <- entrezid
  
  diffex_features <- diffex_features[!is.na(diffex_features$entrezid),]
  rownames(diffex_features) <- diffex_features$entrezid
  
  diffex_features <- new("ExpressionSet", exprs=as.matrix(diffex_features))
  
  # filter by posession of a valid GO term
  haveGO <- sapply(mget(featureNames(diffex_features), org.Hs.egGO), function(x) { if (length(x) == 1 && is.na(x)) FALSE else TRUE})
  diffex_features <- diffex_features[haveGO, ]
  
  ## Format the results for using with GOstats
  diffex_features <- unique(unlist(featureNames(diffex_features)))
  
  return(diffex_features)
  
}

chipEntrezUniverse <- prep_gostats_data(diffex)
q_selectedEntrezIds <- prep_gostats_data(q_data)
# up_q_selectedEntrezIds <- prep_go_data(up_q)
# up_q_selectedEntrezIds <- prep_go_data(dn_q)

chipEntrezUniverse <- prep_topgo_data(diffex)



##
## Run GOstats. I include the conditional, but in practise its not worth much.
##

hgCutoff <- 0.5
pbung <- 0.5



# create topgo object -----------------------------------------------------

sampleGOdata <- new("topGOdata",
                    + description = "Simple session", ontology = "BP",
                    + allGenes = chipEntrezUniverse, geneSel = q_selectedEntrezIds,
                    + nodeSize = 10,
                    + annot = annFUN.db, affyLib = affyLib)

## q selected up/down versus all over 6.25 OVER

params1 <- new("GOHyperGParams", geneIds=q_selectedEntrezIds, universeGeneIds=chipEntrezUniverse, annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
paramsCond1 <- params1
conditional(paramsCond1) <- TRUE
hgOver1 <- hyperGTest(params1)
hgCondOver1 <- hyperGTest(paramsCond1)

## Generate hmtl result files
htmlReport(hgOver1, file=paste("/home/skevin/", dt,"_hgOver_BP_q_v_all_over.html", sep=""))
htmlReport(hgCondOver1, file=paste("/home/skevin/", dt,"_hgCondOver_BP_q_v_all_over.html", sep=""))

# Get the tables out
dat <- as.data.frame(summary(hgOver1))
dat1 <- as.data.frame(summary(hgCondOver1))

write.table(dat, file=paste(dt, "_hgOver_BP_q_v_all_over.txt", sep=""), sep="\t", col.names=T, row.names=T)
write.table(dat1, file=paste(dt, "_hgCondOver_BP_q_v_all_over.txt", sep=""), sep="\t", col.names=T, row.names=T)

# Ok get the graph out

try({

y1 = termGraphs(hgOver1, use.terms=FALSE, pvalue=pbung)
y1.1 = termGraphs(hgCondOver1, use.terms=FALSE, pvalue=pbung)

a1 <- y1$`1`
a2 <- y1.1$`1`

tiff(paste(dt,"_GO_Tree_hgOver_BP_q_v_all_over_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(a1, hgOver1, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

tiff(paste(dt,"_GO_Tree_hgCondOver_BP_q_v_all_over_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(a2, hgCondOver1, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

})

## q selected up/down versus all over 6.25 UNDER

params2 <- new("GOHyperGParams", geneIds=q_selectedEntrezIds, universeGeneIds=chipEntrezUniverse, annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="under")
paramsCond2 <- params2
conditional(paramsCond2) <- TRUE
hgOver2 <- hyperGTest(params2)
hgCondOver2 <- hyperGTest(paramsCond2)

## Generate hmtl result files
htmlReport(hgOver2, file=paste(dt,"_hgOver_BP_q_v_all_under.html", sep=""))
htmlReport(hgCondOver2, file=paste(dt,"_hgCondOver_BP_q_v_all_under.html", sep=""))

# Get the tables out
dat <- as.data.frame(summary(hgOver2))
dat1 <- as.data.frame(summary(hgCondOver2))

write.table(dat, file=paste(dt, "_hgOver_BP_q_v_all_under.txt", sep=""), sep="\t", col.names=T, row.names=T)
write.table(dat1, file=paste(dt, "_hgCondOver_BP_q_v_all_under.txt", sep=""), sep="\t", col.names=T, row.names=T)

# Ok get the graph out

try({

y2 = termGraphs(hgOver2, use.terms=FALSE, pvalue=pbung)
y2.1 = termGraphs(hgCondOver2, use.terms=FALSE, pvalue=pbung)

b1 <- y2$`1`
b2 <- y2.1$`1`

tiff(paste(dt,"_GO_Tree_hgOver_BP_q_v_all_under_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(b1, hgOver2, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

tiff(paste(dt,"_GO_Tree_hgCondOver_BP_q_v_all_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(b2, hgCondOver2, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

})

## q selected up versus all over 6.25 OVER

params3 <- new("GOHyperGParams", geneIds=up_q_selectedEntrezIds, universeGeneIds=chipEntrezUniverse, annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
paramsCond3 <- params3
conditional(paramsCond3) <- TRUE
hgOver3 <- hyperGTest(params3)
hgCondOver3 <- hyperGTest(paramsCond3)

## Generate hmtl result files
htmlReport(hgOver3, file=paste(dt,"_hgOver_BP_up_q_v_all_over.html", sep=""))
htmlReport(hgCondOver3, file=paste(dt,"_hgCondOver_BP_up_q_v_all_over.html", sep=""))

# Get the tables out
dat <- as.data.frame(summary(hgOver3))
dat1 <- as.data.frame(summary(hgCondOver3))

write.table(dat, file=paste(dt, "_hgOver_BP_up_q_v_all_over.txt", sep=""), sep="\t", col.names=T, row.names=T)
write.table(dat1, file=paste(dt, "_hgCondOver_BP_up_q_v_all_over.txt", sep=""), sep="\t", col.names=T, row.names=T)

# Ok get the graph out

try({

y3 = termGraphs(hgOver3, use.terms=FALSE, pvalue=pbung)
y3.1 = termGraphs(hgCondOver3, use.terms=FALSE, pvalue=pbung)

c1 <- y3$`1`
c2 <- y3.1$`1`

tiff(paste(dt,"_GO_Tree_hgOver_BP_up_q_v_all_over_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(c1, hgOver3, node.colors=c(sig="dodgerblue", not="mediumspringgreem"), node.shape="ellipse", add.counts=TRUE)
dev.off()

tiff(paste(dt,"_GO_Tree_hgCondOver_BP_up_q_v_all_over_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(c2, hgCondOver3, node.colors=c(sig="dodgerblue", not="mediumspringgreem"), node.shape="ellipse", add.counts=TRUE)
dev.off()

})

## q selected dn versus all over 6.25 OVER

params4 <- new("GOHyperGParams", geneIds=dn_q_selectedEntrezIds, universeGeneIds=chipEntrezUniverse, annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="over")
paramsCond4 <- params4
conditional(paramsCond4) <- TRUE
hgOver4 <- hyperGTest(params4)
hgCondOver4 <- hyperGTest(paramsCond4)

## Generate hmtl result files
htmlReport(hgOver4, file=paste(dt,"_hgOver_BP_dn_q_v_all_over.html", sep=""))
htmlReport(hgCondOver4, file=paste(dt,"_hgCondOver_BP_dn_q_v_all_over.html", sep=""))

# Get the tables out
dat <- as.data.frame(summary(hgOver4))
dat1 <- as.data.frame(summary(hgCondOver4))

write.table(dat, file=paste(dt, "_hgOver_BP_dn_q_v_all_over.txt", sep=""), sep="\t", col.names=T, row.names=T)
write.table(dat1, file=paste(dt, "_hgCondOver_BP_dn_q_v_all_over.txt", sep=""), sep="\t", col.names=T, row.names=T)

# Ok get the graph out

try({

y4 = termGraphs(hgOver4, use.terms=FALSE, pvalue=pbung)
y4.1 = termGraphs(hgCondOver4, use.terms=FALSE, pvalue=pbung)

d1 <- y4$`1`
d2 <- y4.1$`1`

tiff(paste(dt,"_GO_Tree_hgOver_BP_dn_q_v_all_over_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(d1, hgOver4, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

tiff(paste(dt,"_GO_Tree_hgCondOver_BP_dn_q_v_all_over_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(d2, hgCondOver4, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

})

## q selected up versus all over 6.25 UNDER

params5 <- new("GOHyperGParams", geneIds=up_q_selectedEntrezIds, universeGeneIds=chipEntrezUniverse, annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="under")
paramsCond5 <- params5
conditional(paramsCond5) <- TRUE
hgOver5 <- hyperGTest(params5)
hgCondOver5 <- hyperGTest(paramsCond5)

## Generate hmtl result files
htmlReport(hgOver5, file=paste(dt,"_hgOver_BP_up_q_v_all_under.html", sep=""))
htmlReport(hgCondOver5, file=paste(dt,"_hgCondOver_BP_up_q_v_all_under.html", sep=""))

# Get the tables out
dat <- as.data.frame(summary(hgOver5))
dat1 <- as.data.frame(summary(hgCondOver5))

write.table(dat, file=paste(dt, "_hgOver_BP_up_q_v_all_under.txt", sep=""), sep="\t", col.names=T, row.names=T)
write.table(dat1, file=paste(dt, "_hgCondOver_BP_up_q_v_all_under.txt", sep=""), sep="\t", col.names=T, row.names=T)

# Ok get the graph out

try({

y5 = termGraphs(hgOver5, use.terms=FALSE, pvalue=pbung)
y5.1 = termGraphs(hgCondOver5, use.terms=FALSE, pvalue=pbung)

e1 <- y5$`1`
e2 <- y5.1$`1`

tiff(paste(dt,"_GO_Tree_hgOver_BP_up_q_v_all_under_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(e1, hgOver5, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

tiff(paste(dt,"_GO_Tree_hgCondOver_BP_up_q_v_all_under_0.01.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(e2, hgCondOver5, node.colors=c(sig="dodgerblue", not="mediumspringgreen"), node.shape="ellipse", add.counts=TRUE)
dev.off()

})

## q selected dn versus all over 6.25 UNDER

params6 <- new("GOHyperGParams", geneIds=dn_q_selectedEntrezIds, universeGeneIds=chipEntrezUniverse, annotation="org.Mm.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=FALSE, testDirection="under")
paramsCond6 <- params6
conditional(paramsCond6) <- TRUE
hgOver6 <- hyperGTest(params6)
hgCondOver6 <- hyperGTest(paramsCond6)

## Generate hmtl result files
htmlReport(hgOver6, file=paste(dt,"_hgOver_BP_dn_q_v_all_under.html", sep=""))
htmlReport(hgCondOver6, file=paste(dt,"_hgCondOver_BP_dn_q_v_all_under.html", sep=""))

# Get the tables out
dat <- as.data.frame(summary(hgOver6))
dat1 <- as.data.frame(summary(hgCondOver6))

write.table(dat, file=paste(dt, "_hgOver_BP_dn_q_v_all_under.txt", sep=""), sep="\t", col.names=T, row.names=T)
write.table(dat1, file=paste(dt, "_hgCondOver_BP_dn_q_v_all_under.txt", sep=""), sep="\t", col.names=T, row.names=T)

# Ok get the graph out

try({

y6 = termGraphs(hgOver6, use.terms=FALSE, pvalue=pbung)
y6.1 = termGraphs(hgCondOver6, use.terms=FALSE, pvalue=pbung)

f1 <- y6$`1`
f2 <- y6.1$`1`

tiff(paste(dt,"_GO_Tree_hgOver_BP_dn_q_v_all_under_0.05.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(f1, hgOver6, node.colors=c(sig="dodgerblue", not="mediumspringgreem"), node.shape="ellipse", add.counts=TRUE)
dev.off()

tiff(paste(dt,"_GO_Tree_hgCondOver_BP_dn_q_v_all_under_0.05.tif", sep=""), width=6*1000, height=4*1000)
par(mar=c(5,5,2,2), xaxs ="i", yaxs="i", cex.axis=1.3, cex.lab=1.4)
plotGOTermGraph(f2, hgCondOver6, node.colors=c(sig="dodgerblue", not="mediumspringgreem"), node.shape="ellipse", add.counts=TRUE)
dev.off()

})

save.image(file=paste(dt,"_2.RData", sep=""))

# clean all
rm(list=ls(all=TRUE))

quit("yes")
