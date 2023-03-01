

# load required libraries -------------------------------------------------

library(monocle)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(glue)

# load data
tpm_dat_paths <- c("~/single_cell_projects/quicklinks/FACS_20170407_dshayler_H_sapiens_proj/output/stringtie_transcripts.tpm.csv",
									 "~/single_cell_projects/quicklinks/FACS_20171031_dshayler_H_sapiens_proj/output/stringtie_transcripts.tpm.csv",
									 "~/single_cell_projects/quicklinks/FACS_20181001_dshayler_Organoid_proj/output/transcript_tpm_matrix_fetal.csv",
									 "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/Sunhye_stringtie.tpm.csv",
									 "~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/transcripts.tpm.csv")

# tpm_dat_paths <- c("~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/Sunhye_stringtie.tpm.csv",
# 									"~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/transcripts.tpm.csv")

# readr::write_tsv(tpm_dat, tpm_dat_path)
tpm_dats <- purrr::map(tpm_dat_paths, readr::read_csv)

# load metadata
tpm_meta_paths <- c("~/single_cell_projects/quicklinks/FACS_20170407_dshayler_H_sapiens_proj/Dominic_Cell_Age_and_Sort_Method_042017.csv",
										"~/single_cell_projects/quicklinks/FACS_20171031_dshayler_H_sapiens_proj/FACS_20171031_dshayler_sample_sheet.csv",
										"~/single_cell_projects/quicklinks/FACS_20181001_dshayler_Organoid_proj/10_2018_Seq_3_fetal_cell_metadata.csv",
										"~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/Sunhye_cell_division_day_treatment.csv",
										"~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/FACS_20171031_sunlee_sample_sheet.csv")


tpm_metas <- purrr::map(tpm_meta_paths, readr::read_csv)

expids <- c("ds20170407",
						"ds20171031",
						"ds20181001",
						"shl20170407",
						"shl20171031")

# load useful functions -------------------------------------------------

append_expids <- function(tpm_dat, tpm_meta, expid){
	# browser()
	tpm_ids <- colnames(tpm_dat)[-1]
	meta_ids <- tpm_meta[[1]]
	
	tpm_ids <- str_replace_all(tpm_ids, "[[:alpha:]]", "")
	meta_ids <- str_replace_all(meta_ids, "[[:alpha:]]", "")
	
	tpm_ids <- c("TRANSCRIPT_ID", glue::glue("{expid}_{tpm_ids}"))
	meta_ids <- glue::glue("{expid}_{meta_ids}")
	
	colnames(tpm_dat) <- tpm_ids
	tpm_meta[,1] <- meta_ids
	
	return(list(tpm_dat, tpm_meta))
}

sce_from_tibbles <- function(counts, colData, metadata){
	# browser()
	featuredata <- data.frame(counts[,1])
	rownames(featuredata) <- featuredata[,1]
	
	counts <- data.frame(counts)
	rownames(counts) <- counts[,1]
	counts[,1] <- NULL
	counts <- as.matrix(counts)
	
	colData <- data.frame(colData)
	rownames(colData) <- colData[,1]
	colData <- colData[colnames(counts),]
	sumexp <- SummarizedExperiment(assays=list(counts=counts), colData=colData, rowData=featuredata, metadata=metadata)	
	
	sce <- as(sumexp, "SingleCellExperiment")
	return(sce)
}



mergeSingleCellExperiments <- function(sce1, sce2){
	# browser()
	
	# check that features (genes or transcripts) are idential between sces
	identical(rownames(sce1), rownames(sce2))
	
	#	Now we’ll check that there aren’t any repeated cellIDs:
	
	sum(colnames(sce1) %in% colnames(sce2))
	
	#	Everything is ok, so we can go ahead and combine them:
	
	sample_ids <- c(colnames(sce1), colnames(sce2))
	merge_cd <- dplyr::bind_rows(data.frame(colData(sce1)), data.frame(colData(sce2)))
	rownames(merge_cd) <- sample_ids
	if(!"Sample_ID" %in% colnames(merge_cd)){
		merge_cd <- cbind(Sample_ID = rownames(merge_cd), merge_cd)
	}
	
	# merge counts
	merge_counts <- merge(counts(sce1), counts(sce2), by = 0)
	rownames(merge_counts) <- merge_counts[,1]
	merge_counts[,1] <- NULL
	
	merge_rd <- data.frame(row.names = rownames(merge_counts), "features" = rownames(merge_counts))
	
	merge_counts <- as.matrix(merge_counts)	
	
	merge_sce <- SingleCellExperiment(assays=list(counts=merge_counts), colData=merge_cd, rowData=merge_rd)
	
	return(merge_sce)
}



# make singlecellexperiments ----------------------------------------------

exp_lists <- purrr::pmap(list(tpm_dats, tpm_metas, expids), append_expids)

counts <- lapply(exp_lists, "[[", 1)
metas <- lapply(exp_lists, "[[", 2)

sces <- purrr::map2(counts, metas, sce_from_tibbles, list(c("test")))

merge_sce <- purrr::reduce(sces, mergeSingleCellExperiments)

run_census <- function(sce, census_output_file){
		
	# create new celldataset
	pd <- new("AnnotatedDataFrame", data=data.frame(colData(sce)))
	fd <- new("AnnotatedDataFrame", data=data.frame(rowData(sce)))
	
	
	HSMM <- newCellDataSet(counts(sce),
												 phenoData = pd,
												 featureData = fd,
												 lowerDetectionLimit=0.1,
												 expressionFamily=tobit(Lower=0.1))

	HSMM_path <- gsub("census_matrix.csv", "raw_dataset.rds", census_output_file)
	
	saveRDS(HSMM, HSMM_path)
	
	rpc_matrix <- monocle::relative2abs(HSMM, method = "num_genes")
	
	
	# Now, make a new CellDataSet using the RNA counts
	HSMM <- newCellDataSet(as(rpc_matrix, "sparseMatrix"),
												 phenoData = pd,
												 featureData = fd,
												 lowerDetectionLimit=1,
												 expressionFamily=negbinomial.size())
	
	
	# HSMM <- estimateSizeFactors(HSMM)
	# HSMM <- estimateDispersions(HSMM)
	
	
	#print output of census to csv prior to monocle workflow
	
	write.csv(as.matrix(Biobase::exprs(HSMM)), census_output_file)
	
	census_meta_file <- gsub("census_matrix.csv", "census_meta.csv", census_output_file)

	write.csv(pData(HSMM), census_meta_file, row.names = FALSE)
	
	return(as.matrix(Biobase::exprs(HSMM)))

}

census_path <- "~/single_cell_projects/quicklinks/5_seq_proj/output/transcripts_census_matrix.csv"
census_mat <- run_census(merge_sce, census_path)
