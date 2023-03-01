#!/usr/bin/Rscript

library(recount)

## Find a project of interest
project_info <- abstract_search('rod cell')

## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)

## Load the data
load(file.path(project_info$project, 'rse_gene.Rdata'))

## Browse the project at SRA
browse_study(project_info$project)

## View GEO ids
colData(rse_gene)$geo_accession

## Extract the sample characteristics
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)


## We can now define some sample information to use
sample_info <- data.frame(
  run = colData(rse_gene)$run,
  group = ifelse(grepl('GFP negative', colData(rse_gene)$characteristics), 'positive', 'induced'),
  gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
                                                                                'targeting ')[[1]][2], ' gene')[[1]][1] }),
  cell.line = geochar$cell.line
)