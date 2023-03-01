#!/usr/bin/Rscript


## ---- echo=FALSE, results='hide'-------------------------------------------
all.urls <- c("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81076/suppl/GSE81076%5FD2%5F3%5F7%5F10%5F17%2Etxt%2Egz", 
              "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz", 
              "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/files/E-MTAB-5061.processed.1.zip",
              "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt")

all.basenames <- basename(all.urls)
all.basenames[1] <- "GSE81076_D2_3_7_10_17.txt.gz" 
all.basenames[2] <- "GSE85241_cellsystems_dataset_4donors_updated.csv"

all.modes <- c("wb", "w", "wb", "w")
for (x in seq_along(all.urls)) { 
  if (!file.exists(all.basenames[x])) {
    download.file(all.urls[x], all.basenames[x], mode=all.modes[x])
  }
}

## --------------------------------------------------------------------------
gse81076.df <- read.table("GSE81076_D2_3_7_10_17.txt.gz", sep='\t', 
                          header=TRUE, stringsAsFactors=FALSE, row.names=1)
dim(gse81076.df)

## --------------------------------------------------------------------------
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse81076.df))
table(donor.names)
plate.id <- sub("^D[0-9]+(.*)_.*", "\\1", colnames(gse81076.df))
table(plate.id)

## --------------------------------------------------------------------------
gene.symb <- gsub("__chr.*$", "", rownames(gse81076.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)

library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.symb, keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.symb[is.spike]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
gse81076.df <- gse81076.df[keep,]
rownames(gse81076.df) <- gene.ids[keep]
summary(keep)

## --------------------------------------------------------------------------
library(SingleCellExperiment)
sce.gse81076 <- SingleCellExperiment(list(counts=as.matrix(gse81076.df)),
                                     colData=DataFrame(Donor=donor.names, Plate=plate.id),
                                     rowData=DataFrame(Symbol=gene.symb[keep]))
isSpike(sce.gse81076, "ERCC") <- grepl("^ERCC-", rownames(gse81076.df)) 
sce.gse81076  

## --------------------------------------------------------------------------
library(scater)
sce.gse81076 <- calculateQCMetrics(sce.gse81076, compact=TRUE)
QC <- sce.gse81076$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE))

## --------------------------------------------------------------------------
discard <- low.lib | low.genes | high.spike
sce.gse81076 <- sce.gse81076[,!discard]
summary(discard)

## --------------------------------------------------------------------------
library(scran)
clusters <- quickCluster(sce.gse81076, min.mean=0.1)
table(clusters)
sce.gse81076 <- computeSumFactors(sce.gse81076, min.mean=0.1, clusters=clusters)
summary(sizeFactors(sce.gse81076))

## --------------------------------------------------------------------------
sce.gse81076 <- computeSpikeFactors(sce.gse81076, general.use=FALSE)
summary(sizeFactors(sce.gse81076, "ERCC"))

## --------------------------------------------------------------------------
sce.gse81076 <- normalize(sce.gse81076)

## ----var-gse81076, fig.cap="Variance of normalized log-expression values for each gene in the GSE81076 dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
block <- paste0(sce.gse81076$Plate, "_", sce.gse81076$Donor)
fit <- trendVar(sce.gse81076, block=block, parametric=TRUE) 
dec <- decomposeVar(sce.gse81076, fit)
plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
OBis.spike <- isSpike(sce.gse81076)
points(dec$mean[is.spike], dec$total[is.spike], col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

## --------------------------------------------------------------------------
dec.gse81076 <- dec
dec.gse81076$Symbol <- rowData(sce.gse81076)$Symbol
dec.gse81076 <- dec.gse81076[order(dec.gse81076$bio, decreasing=TRUE),]
head(dec.gse81076)

## ---- echo=FALSE, results="hide"-------------------------------------------
rm(gse81076.df)
gc()

## --------------------------------------------------------------------------
gse85241.df <- read.table("GSE85241_cellsystems_dataset_4donors_updated.csv", 
                          sep='\t', h=TRUE, row.names=1, stringsAsFactors=FALSE)
dim(gse85241.df)

## --------------------------------------------------------------------------
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse85241.df))
table(donor.names)
plate.id <- sub("^D[0-9]+\\.([0-9]+)_.*", "\\1", colnames(gse85241.df))
table(plate.id)

## --------------------------------------------------------------------------
gene.symb <- gsub("__chr.*$", "", rownames(gse85241.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)

library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.symb, keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.symb[is.spike]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
gse85241.df <- gse85241.df[keep,]
rownames(gse85241.df) <- gene.ids[keep]
summary(keep)

## --------------------------------------------------------------------------
sce.gse85241 <- SingleCellExperiment(list(counts=as.matrix(gse85241.df)),
                                     colData=DataFrame(Donor=donor.names, Plate=plate.id),
                                     rowData=DataFrame(Symbol=gene.symb[keep]))
isSpike(sce.gse85241, "ERCC") <- grepl("^ERCC-", rownames(gse85241.df)) 
sce.gse85241  

## --------------------------------------------------------------------------
sce.gse85241 <- calculateQCMetrics(sce.gse85241, compact=TRUE)
QC <- sce.gse85241$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE))

## --------------------------------------------------------------------------
discard <- low.lib | low.genes | high.spike
sce.gse85241 <- sce.gse85241[,!discard]
summary(discard)

## --------------------------------------------------------------------------
clusters <- quickCluster(sce.gse85241, min.mean=0.1, method="igraph")
table(clusters)
sce.gse85241 <- computeSumFactors(sce.gse85241, min.mean=0.1, clusters=clusters)
summary(sizeFactors(sce.gse85241))
sce.gse85241 <- computeSpikeFactors(sce.gse85241, general.use=FALSE)
summary(sizeFactors(sce.gse85241, "ERCC"))
sce.gse85241 <- normalize(sce.gse85241)

## ----var-gse85241, fig.cap="Variance of normalized log-expression values for each gene in the GSE85241 dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red)."----
block <- paste0(sce.gse85241$Plate, "_", sce.gse85241$Donor)
fit <- trendVar(sce.gse85241, block=block, parametric=TRUE) 
dec <- decomposeVar(sce.gse85241, fit)
plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
is.spike <- isSpike(sce.gse85241)
points(dec$mean[is.spike], dec$total[is.spike], col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

## --------------------------------------------------------------------------
dec.gse85241 <- dec
dec.gse85241$Symbol <- rowData(sce.gse85241)$Symbol
dec.gse85241 <- dec.gse85241[order(dec.gse85241$bio, decreasing=TRUE),]
head(dec.gse85241)

## ---- echo=FALSE, results="hide"-------------------------------------------
rm(gse85241.df)
gc()

## --------------------------------------------------------------------------
unzip("E-MTAB-5061.processed.1.zip")

# Figuring out the number of libraries (-1 for the '#sample').
header <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt", 
                     nrow=1, sep="\t", comment.char="", stringsAsFactors=FALSE)
ncells <- ncol(header) - 1L

# Loading only the gene names and the counts.
col.types <- vector("list", ncells*2 + 2)
col.types[1:2] <- "character"
col.types[2+ncells + seq_len(ncells)] <- "integer"
e5601.df <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt", 
                       sep="\t", colClasses=col.types)

# Disentangling the gene names and the counts.
gene.data <- e5601.df[,1:2]
e5601.df <- e5601.df[,-(1:2)]
colnames(e5601.df) <- as.character(header[1,-1])
dim(e5601.df)

## --------------------------------------------------------------------------
is.spike <- grepl("^ERCC-", gene.data[,2])
table(is.spike)

library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.data[,1], keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.data[is.spike,2]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
e5601.df <- e5601.df[keep,]
rownames(e5601.df) <- gene.ids[keep]
summary(keep)

## --------------------------------------------------------------------------
metadata <- read.table("E-MTAB-5061.sdrf.txt", header=TRUE, 
                       sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
m <- match(colnames(e5601.df), metadata[["Assay Name"]])
stopifnot(all(!is.na(m)))
metadata <- metadata[m,]
donor.id <- metadata[["Characteristics[individual]"]]
table(donor.id)

## --------------------------------------------------------------------------
sce.e5601 <- SingleCellExperiment(list(counts=as.matrix(e5601.df)),
                                  colData=DataFrame(Donor=donor.id),
                                  rowData=DataFrame(Symbol=gene.data[keep,1]))
isSpike(sce.e5601, "ERCC") <- grepl("^ERCC-", rownames(e5601.df)) 
sce.e5601  

## --------------------------------------------------------------------------
sce.e5601 <- calculateQCMetrics(sce.e5601, compact=TRUE)
QC <- sce.e5601$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3) 
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
low.spike <- isOutlier(QC$feature_control_ERCC$log10_total_counts, type="lower", nmad=2)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE), LowSpike=sum(low.spike))

## --------------------------------------------------------------------------
discard <- low.lib | low.genes | high.spike | low.spike
sce.e5601 <- sce.e5601[,!discard]
summary(discard)

## --------------------------------------------------------------------------
clusters <- quickCluster(sce.e5601, min.mean=1, method="igraph")
table(clusters)
sce.e5601 <- computeSumFactors(sce.e5601, min.mean=1, clusters=clusters)
summary(sizeFactors(sce.e5601))
sce.e5601 <- computeSpikeFactors(sce.e5601, general.use=FALSE)
summary(sizeFactors(sce.e5601, "ERCC"))
sce.e5601 <- normalize(sce.e5601)

## ----var-e5601, fig.cap="Variance of normalized log-expression values for each gene in the E-MTAB-5601 dataset, plotted against the mean log-expression. Each plot corresponds to a donor, where the blue line represents the mean-dependent trend fitted to the variances of the spike-in transcripts (red).", fig.width=6, fig.asp=2.5----
donors <- sort(unique(sce.e5601$Donor))
is.spike <- isSpike(sce.e5601)
par(mfrow=c(ceiling(length(donors)/2), 2), 
    mar=c(4.1, 4.1, 2.1, 0.1))
collected <- list()
for (x in unique(sce.e5601$Donor)) {
  current <- sce.e5601[,sce.e5601$Donor==x]
  if (ncol(current)<2L) { next }
  current <- normalize(current)
  fit <- trendVar(current, parametric=TRUE) 
  dec <- decomposeVar(current, fit)
  plot(dec$mean, dec$total, xlab="Mean log-expression",
       ylab="Variance of log-expression", pch=16, main=x)
  points(fit$mean, fit$var, col="red", pch=16)
  curve(fit$trend(x), col="dodgerblue", add=TRUE)
  collected[[x]] <- dec
}

## --------------------------------------------------------------------------
dec.e5601 <- do.call(combineVar, collected)
dec.e5601$Symbol <- rowData(sce.e5601)$Symbol
dec.e5601 <- dec.e5601[order(dec.e5601$bio, decreasing=TRUE),]
head(dec.e5601)

## ---- echo=FALSE, results="hide"-------------------------------------------
rm(e5601.df)
gc()

## --------------------------------------------------------------------------
top.e5601 <- rownames(dec.e5601)[seq_len(1000)]
top.gse85241 <- rownames(dec.gse85241)[seq_len(1000)]
top.gse81076 <- rownames(dec.gse81076)[seq_len(1000)]
chosen <- Reduce(intersect, list(top.e5601, top.gse85241, top.gse81076))

# Adding some gene symbols for interpretation.
symb <- mapIds(org.Hs.eg.db, keys=chosen, keytype="ENSEMBL", column="SYMBOL")
DataFrame(ID=chosen, Symbol=symb)

## --------------------------------------------------------------------------
# Identifying genes that are annotated in all batches.
in.all <- Reduce(intersect, list(rownames(dec.e5601), 
                                 rownames(dec.gse85241), rownames(dec.gse81076)))

# Setting weighted=FALSE so each batch contributes equally.
combined <- combineVar(dec.e5601[in.all,], dec.gse85241[in.all,],
                       dec.gse81076[in.all,], weighted=FALSE)
chosen2 <- rownames(combined)[head(order(combined$bio, decreasing=TRUE), 1000)]

## --------------------------------------------------------------------------
original <- list(logcounts(sce.e5601)[chosen,],
                 logcounts(sce.gse81076)[chosen,],
                 logcounts(sce.gse85241)[chosen,])
corrected <- do.call(mnnCorrect, c(original, list(k=20, sigma=0.1)))
str(corrected$corrected)

## --------------------------------------------------------------------------
corrected$pairs

## --------------------------------------------------------------------------
omat <- do.call(cbind, original)
mat <- do.call(cbind, corrected$corrected)
colnames(mat) <- NULL
sce <- SingleCellExperiment(list(original=omat, corrected=mat))
colData(sce)$Batch <- rep(c("e5601", "gse81076", "gse85241"),
                          lapply(corrected$corrected, ncol))
sce

## ----tsne-batch, fig.width=10, fig.asp=0.6, fig.cap="t-SNE plots of the pancreas datasets, before and after MNN correction. Each point represents a cell and is coloured by the batch of origin."----
osce <- runTSNE(sce, exprs_values="original", rand_seed=100)
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")
csce <- runTSNE(sce, exprs_values="corrected", rand_seed=100)
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
multiplot(ot, ct, cols=2)

## ----tsne-markers, fig.width=10, fig.height=10, fig.cap="t-SNE plots after MNN correction, where each point represents a cell and is coloured by its corrected expression of key marker genes for known cell types in the pancreas."----
ct.gcg <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000115263") + ggtitle("Alpha cells (GCG)")
ct.ins <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000254647") + ggtitle("Beta cells (INS)")
ct.sst <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000157005") + ggtitle("Delta cells (SST)")
ct.ppy <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000108849") + ggtitle("PP cells (PPY)")
multiplot(ct.gcg, ct.ins, ct.sst, ct.ppy, cols=2)
