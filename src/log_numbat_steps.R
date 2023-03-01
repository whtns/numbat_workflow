#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

resources_path <- "/dataVolume/storage/single_cell_projects/resources/"

combined_metadata <- read_tsv("data/combined_metadata.tsv") %>%
  dplyr::mutate(numbat_log = fs::path(resources_path, paste0(study, "_et_al_proj"), "output/numbat/", run, "log.txt")) %>%
  dplyr::mutate(numbat_log_exists = file_exists(numbat_log))

read_numbat_log <- function(numbat_log_file){
  browser()
  test0 <- read_lines(numbat_log_file) %>%
    str_subset(" = ") %>%
    tibble::as_tibble() %>%
    tidyr::separate(value, c("parameter", "value"), sep = " = ")
}

test0 <-
  combined_metadata %>%
  dplyr::filter(file_exists(numbat_log)) %>%
  dplyr::pull(numbat_log) %>%
  purrr::map(read_numbat_log) %>%
  identity()

test0 <- purrr::map(combined_metadata$numbat_log, ~possibly(read_numbat_log(.x)))

numbat_logs
