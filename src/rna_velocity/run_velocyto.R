#!/usr/bin/env Rscript

## load packags ------------------------------------------------------------------------
library(tidyverse)
library(velocyto.R)
library(Seurat)
library(rprojroot)
library(fs)
library(gsubfn)
proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))


## load objects------------------------------------------------------------------------

set_loom_names <- function(loom_mats, batch){
	# browser()
  oldnames <- colnames(loom_mats[[1]])
  
  toreplace<-list("MyTissue:" = paste0(batch, "_S"), "_removed_duplicates.bam" = "")
  newnames <- gsubfn(paste(names(toreplace),collapse="|"),toreplace,oldnames)

	name_loom_mat <- function(loom_mat, newnames){
		colnames(loom_mat) <- newnames
		return(loom_mat)
	}
	
	named_loom_mats <- purrr::map(loom_mats, name_loom_mat, newnames)	
	return(named_loom_mats)
}

save_anno <- function(seu, path){
	anno_tbl <- as_tibble(seu[[]])
	write_csv(anno_tbl, path)
}


# subset ldat by seurat object.size
ldat <- purrr::map(velocity_params$ldat, ~.x[,colnames(.x) %in% colnames(velocity_params$seu)])

# subset seurat object by ldat
seu <- velocity_params$seu[,colnames(velocity_params$seu) %in% colnames(ldat[[1]])]

## grab cell colors ------------------------------------------------------------------------

p <- DimPlot(seu, reduction = "umap", group.by = "clusters_1")
col_vec <- unique(ggplot_build(p)$data[[1]]) %>% 
  arrange(group) %>% 
  pull(colour) %>% 
  unique()

DimPlot(seu, reduction = "umap", group.by = "clusters_1")
DimPlot(seu, reduction = "umap", cols = col_vec, group.by = "clusters_1")
# DimPlot(seu, reduction = "umap", cols = col_vec)


## format cell colors------------------------------------------------------------------------

cell.colors <- as_tibble(seu[["clusters_1"]], rownames = "cellid") %>% 
  tibble::deframe() %>% 
  as.factor()

# levels(cell.colors) <- as.character(as.numeric(levels(cell.colors)) + 1)
levels(cell.colors) <- col_vec

# prep embedding ------------------------------
emb <- Embeddings(seu, velocity_params$emb_flag)


## organzie loom data ------------------------------------------------------------------------
# exonic read (spliced) expression matrix
emat <- ldat$spliced;
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning;
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))


## ------------------------------------------------------------------------
# and if we use spanning reads (smat)
length(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))



# plot_spliced_mag <- function(ldat) {
#   hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
# }
# 
# pdf(fs::path(batch_dir, "reads_per_gene.pdf"))
# plot_spliced_mag(ldat)
# dev.off()



# 
# ## calculate velocity------------------------------------------------------------------------
fit.quantile <- 0.05;
rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 5, fit.quantile = fit.quantile)

# save velocity to seurat object in slot @misc$vel
velocity_params$seu@misc$vel <- rvel.qf
saveRDS(velocity_params$seu, velocity_params$seu_path)
