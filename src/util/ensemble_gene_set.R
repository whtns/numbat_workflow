#!/usr/bin/env Rscript


# load libraries ----------------------------------------------------------
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)


# load reference ---------------------------------------------------------------
m_df = msigdbr(species = "Homo sapiens")

m_t2g = m_df %>% 
	dplyr::select(gs_name, entrez_gene) %>% 
	as.data.frame()


# load data ---------------------------------------------------------------

shl_data <- read_csv("~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/scde/transcript/20181025_PT1_scde_shCtrl_vs_shRB1-LatePT1.csv")


# run fgsea ---------------------------------------------------------------

# rank gene list

rank_feature_list <- function(feature_df){
	
	feature_df <- dplyr::mutate(feature_df, fcSign = sign(mle)) %>% 
		dplyr::mutate(rank = -log10(cZ)/fcSign) %>% 
		dplyr::select(symbol, rank)
	
	feature_df <- feature_df[complete.cases(feature_df),]
	
	ranked_list <- feature_df[["rank"]]
	names(ranked_list) <- feature_df[["symbol"]]
	
	return(ranked_list)

}

ranked_list <- rank_feature_list(shl_data)


fgseaRes <- fgsea(pathways = examplePathways, 
									stats = ranked_list,
									minSize=15,
									maxSize=500,
									nperm=10000)

m_list = m_df %>% 
	split(x = .$gene_symbol, f = .$gs_name)

fgsea(pathways = m_list, ...)



clusterProfiler::enricher(gene = genes_entrez, TERM2GENE = m_t2g, ...)
