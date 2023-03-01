#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

collin_metadata <- read_csv("data/external_rb_scrnaseq_metadata/collin_et_al_metadata.csv") %>%
  dplyr::mutate(Run = sample_id)

field_metadata <- read_csv("data/external_rb_scrnaseq_metadata/field_et_al_metadata(1).txt")

wu_metadata <- read_csv("data/external_rb_scrnaseq_metadata/wu_et_al_metadata.csv") %>%
  dplyr::rename(`SRA Study` = SRAStudy)

yang_metadata <- read_csv("data/external_rb_scrnaseq_metadata/yang_et_al_metadata.txt")

combined_colnames <- c("study", "sample_id", "names", "type", "age_mo", "assay_type",
                       "avg_spot_len", "bases", "bio_project", "bio_sample", "bytes",
                       "center_name", "consent", "datastore_filetype", "datastore_provider",
                       "datastore_region", "developmental_stage", "experiment", "geo_accession_exp",
                       "instrument", "library_layout", "library_selection", "library_source",
                       "organism", "platform", "release_date", "sample_name", "source_name",
                       "sra_study", "tissue", "run", "familial", "iirc_tumor_class",
                       "library_name", "load_date", "spots", "spots_with_mates", "avg_length",
                       "size_mb", "assembly_name", "download_path", "library_strategy",
                       "insert_size", "insert_dev", "model", "study_pubmed_id", "project_id",
                       "sample", "sample_type", "tax_id", "scientific_name", "g1k_pop_code",
                       "source", "g1k_analysis_group", "subject_id", "sex", "disease",
                       "tumor", "affection_status", "analyte_type", "histological_type",
                       "body_site", "submission", "dbgap_study_accession", "run_hash",
                       "read_hash", "age", "biomaterial_provider", "bio_sample_model",
                       "isolate", "replicate", "treatment")

combined_metadata <-
  list(collin = collin_metadata,
       field = field_metadata,
       wu = wu_metadata,
       yang = yang_metadata) %>%
  purrr::map(janitor::clean_names) %>%
  dplyr::bind_rows(.id = "study")
