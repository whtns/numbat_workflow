#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

metadata <- read_csv("data/20211112-SHL-FACS-Hs_metadata.csv") %>%
    # janitor::clean_names() %>%
    # dplyr::select(-sample_id) %>%
    # dplyr::rename(sample_id = sample_name) %>%
    # dplyr::select(-x11) %>% 
    dplyr::mutate(sample_number = str_extract(sample_id, "[0-9]+")) %>%
    dplyr::mutate(sample_number = str_pad(sample_number, width = max(nchar(sample_number)), pad = "0")) %>%
    dplyr::mutate(sample_id = paste0("20211112-SHL-FACS-Hs", "-", sample_number)) %>% 
    dplyr::mutate(names = sample_id, type = "PE") %>%
    identity()

write_tsv(metadata, "data/metadata.txt")

write_csv(metadata, "data/20211112-SHL-FACS-Hs_metadata.csv")

proj_name <- here::here() %>% 
    fs::path_file() %>% 
    str_remove("_proj")

fastq_files <- 
    dir_ls("data/", glob = "*.gz", recurse = TRUE) %>% 
    tibble::enframe("name", "path") %>% 
    dplyr::mutate(pair = dplyr::case_when(str_detect(name, "R1") ~ "R1",
                                          str_detect(name, "R2") ~ "R2",)) %>%
    dplyr::mutate(new_path = paste0(str_remove(name, "_.*"), "_", pair, ".fastq.gz")) %>%
    dplyr::mutate(sample_number = str_extract(new_path, "[0-9]+_")) %>% 
    dplyr::mutate(sample_number = str_remove(sample_number, "_")) %>% 
    dplyr::mutate(sample_number = str_pad(sample_number, width = max(nchar(sample_number)), pad = "0")) %>%
    dplyr::mutate(new_path = paste0(proj_name, "-", sample_number, "_", pair, ".fastq.gz")) %>%
    dplyr::mutate(new_path = fs::path(path_dir(path), new_path)) %>%
    identity()

purrr::map2(fastq_files$path, fastq_files$new_path, fs::file_move)
