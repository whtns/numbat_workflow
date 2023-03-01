#!/usr/bin/Rscript


# load required libraries -------------------------------------------------

library(shiny)
library(heatmaply)
library(shinyHeatmaply)

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


test <- reorder_mat(mat_path, "rod_markers", "NRL")

