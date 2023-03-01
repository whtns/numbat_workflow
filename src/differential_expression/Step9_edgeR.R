#!/usr/bin/Rscript --vanilla --slave

## Change to data directory
# setwd("/data/met/")

today <- Sys.Date()
dt <- format(today, format="%d%b%y")

library(RUVSeq)
library(RColorBrewer)

## import and merge the output of htseq-count

KO1 <- read.csv("sorted_18Feb16_WS_YL04-KO1_gencode_m8_ERCC_star.Aligned.out.txt", sep="", header=F)
colnames(KO1) <- c("EnsemblID", "Set18_KO1")
KO2 <- read.csv("sorted_18Feb16_WS_YL05-KO2_gencode_m8_ERCC_star.Aligned.out.txt", sep="", header=F)
colnames(KO2) <- c("EnsemblID", "Set18_KO2")
KO3 <- read.csv("sorted_18Feb16_WS_YL06-KO3_gencode_m8_ERCC_star.Aligned.out.txt", sep="", header=F)
colnames(KO3) <- c("EnsemblID", "Set18_KO3")
WT1 <- read.csv("sorted_18Feb16_WS_YL01-WT1_gencode_m8_ERCC_star.Aligned.out.txt", sep="", header=F)
colnames(WT1) <- c("EnsemblID", "Set18_WT1")
WT2 <- read.csv("sorted_18Feb16_WS_YL02-WT2_gencode_m8_ERCC_star.Aligned.out.txt", sep="", header=F)
colnames(WT2) <- c("EnsemblID", "Set18_WT2")
WT3 <- read.csv("sorted_18Feb16_WS_YL03-WT3_gencode_m8_ERCC_star.Aligned.out.txt", sep="", header=F)
colnames(WT3) <- c("EnsemblID", "Set18_WT3")
one <- merge(KO1, KO2, by="EnsemblID")
two <- merge(one, KO3, by="EnsemblID")
three <- merge(two, WT1, by="EnsemblID")
four <- merge(three, WT2, by="EnsemblID")
five <- merge(four, WT3, by="EnsemblID")
write.table(five, file=paste(dt, "_WS_Set18_raw_merged.txt", sep=""), sep="\t", col.names=T, row.names=F)

## This is when you should quit and filter. Or put in a better function

dat <- five
len <- read.csv("gencode_M8_lengths.txt", sep="", header=F)
colnames(len) <- c("EnsemblID", "length", "NoExons")
data <- merge(dat, len, by="EnsemblID", all.x=T)
write.table(data, file=paste(dt, "_Set18_raw_lengths.txt", sep=""), sep="\t", col.names=T, row.names=F)

## Get the RPKM values for the data, note calcNormFactors is an essential step

rownames(data) <- data$EnsemblID
keeps <- c("Set18_KO1", "Set18_KO2", "Set18_KO3", "Set18_WT1", "Set18_WT2", "Set18_WT3", "length")
data2 <- data[keeps]

y <- DGEList(counts=data2[,1:6], genes=data.frame(Length=data2[,7]))
y <- calcNormFactors(y)
RPKM <- rpkm(y)
colnames(RPKM) <- c("Set18_KO1_RPKM", "Set18_KO2_RPKM", "Set18_KO3_RPKM", "Set18_WT1_RPKM", "Set18_WT2_RPKM", "Set18_WT3_RPKM")
one <- merge(data, RPKM, by=0)
write.table(one, file=paste(dt, "_Set18_raw_length_rpkm.txt", sep=""), sep="\t", col.names=T, row.names=F)

## Now you should probably quit and make kernel density plots for additional filtering

## Process data using RUVSeq and the ERCC controls. It is a slightly different proceedure if you are using both mixes. In which case you would manually specify only the Group B probes as they are the same concentration between both mixes. In the below case only ERCC mix was used for all six samples (not ideal)

keeps <- c("Set18_KO1", "Set18_KO2", "Set18_KO3", "Set18_WT1", "Set18_WT2", "Set18_WT3")
data4 <- data2[keeps]
spikes <- rownames(data4)[grep("^ERCC", rownames(data4))]
x <- as.factor(c("KO", "KO", "KO", "CTR", "CTR", "CTR"))
set <- newSeqExpressionSet(as.matrix(data4), phenoData = data.frame(x, row.names=colnames(data4)))

colors <- brewer.pal(3, "Set1")

tiff(file=paste(dt, "_RawData_RLE.tiff", sep=""), width=5*300, height=5*300, res=300)
plotRLE(set, outline=FALSE, ylim=c(-4,4), col=colors[x])
dev.off()

tiff(file=paste(dt, "_RawData_RLE.tiff", sep=""), width=5*300, height=5*300, res=300)
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()

set1 <- betweenLaneNormalization(set, which="upper")

tiff(file=paste(dt, "_UQ_Norm_RLE.tiff", sep=""), width=5*300, height=5*300, res=300)
plotRLE(set1, outline=FALSE, ylim=c(-1.5,1.5), col=colors[x])
dev.off()

tiff(file=paste(dt, "_UQ_Norm_PCA.tiff", sep=""), width=5*300, height=5*300, res=300)
plotPCA(set1, col=colors[x], cex=0.8)
dev.off()

set2 <- RUVg(set, spikes, k=1)

tiff(file=paste(dt, "_RUVg_Spikes_RLE.tiff", sep=""), width=5*300, height=5*300, res=300)
plotRLE(set2, outline=FALSE, ylim=c(-1.5,1.5), col=colors[x])
dev.off()

tiff(file=paste(dt, "_RUVg_Spikes_PCA.tiff", sep=""), width=5*300, height=5*300, res=300)
plotPCA(set2, col=colors[x], cex=0.8)
dev.off()

##  Now process the data with edgeR. There are a ton of different design matrices for edgeR same as limma.

## Process with no normalization at all log will contain numbers of de genes and the BCV. The BCV should be minimized
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design, verbose=T)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]

tiff(file=paste(dt, "_BCVgraph_NoNormalization.tiff", sep=""), width=8*300, height=5*300, res=300)
plotBCV(y)
dev.off()

tiff(file=paste(dt, "_Smear_NoNormalization.tiff", sep=""), width=8*300, height=5*300, res=300)
plotSmear(lrt, de.tags=detags)
abline(h=c(-2, 2), col="blue")
dev.off()

proc1 <- merge(top, data, by=0, all.y=T)
write.table(proc1, file=paste(dt, "_ManuallyFiltered_EdgeR_NoNormalization.txt", sep=""), sep="\t", col.names=T, row.names=F)

## Use the un-normalized data for an empirical normalization using non DE genes
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
set3 <- RUVg(set, empirical, k=1)

tiff(file=paste(dt, "_RUVg_Empirical_RLE.tiff", sep=""), width=5*300, height=5*300, res=300)
plotRLE(set3, outline=FALSE, ylim=c(-4,4), col=colors[x])
dev.off()

tiff(file=paste(dt, "_RUVg_Empirical_PCA.tiff", sep=""), width=5*300, height=5*300, res=300)
plotPCA(set3, col=colors[x], cex=1.2)
dev.off()

## Now process with the 'spikes' set and collect BCV and differential gene expression
design <- model.matrix(~x + W_1, pData(set2))
y <- DGEList(counts=normCounts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design, verbose=T)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top1 <- topTags(lrt, n=nrow(set))$table
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]

tiff(file=paste(dt, "_BCVgraph_Spikes.tiff", sep=""), width=8*300, height=5*300, res=300)
plotBCV(y)
dev.off()

tiff(file=paste(dt, "_Smear_Spikes.tiff", sep=""), width=8*300, height=5*300, res=300)
plotSmear(lrt, de.tags=detags)
abline(h=c(-2, 2), col="blue")
dev.off()

proc2 <- merge(top1, data, by=0, all.y=T)
write.table(proc1, file=paste(dt, "_ManuallyFiltered_EdgeR_Spikes.txt", sep=""), sep="\t", col.names=T, row.names=F)

## Now process with the 'Empirical' data set collect BCV and DE genes
design <- model.matrix(~x + W_1, pData(set3))
y <- DGEList(counts=normCounts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design, verbose=T)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top1 <- topTags(lrt, n=nrow(set))$table
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]

tiff(file=paste(dt, "_BCVgraph_Empirical.tiff", sep=""), width=8*300, height=5*300, res=300)
plotBCV(y)
dev.off()

tiff(file=paste(dt, "_Smear_Empirical.tiff", sep=""), width=8*300, height=5*300, res=300)
plotSmear(lrt, de.tags=detags)
abline(h=c(-2, 2), col="blue")
dev.off()

proc2 <- merge(top1, data, by=0, all.y=T)
write.table(proc1, file=paste(dt, "_ManuallyFiltered_EdgeR_Empirical.txt", sep=""), sep="\t", col.names=T, row.names=F)

# Save an image
save.image(file=paste(dt,"_1.RData". sep=""))

# clean all
rm(list=ls(all=TRUE))

quit("yes")

