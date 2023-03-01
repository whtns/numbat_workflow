#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(glue)
library(numbat)

source("functions.R")

numbat_dirs <- dir_ls("output/numbat/", type = "directory")

safe_make_numbat_plots <- purrr::safely(make_all_numbat_plots)

map(numbat_dirs, safe_make_numbat_plots)

# safe_make_numbat_plots("output/numbat/SRR13884240")
