fetal_cell_path <- "10_2018_Seq_3_Fetal_cell_metadata.csv"
fetal_cell_meta2 <- read.csv(fetal_cell_path, row.names = 1) %>% 
	dplyr::mutate(Sample_ID = gsub("_", "", Sample_ID))
write.csv(fetal_cell_meta, fetal_cell_path)


all_cell_path <- "10_2018_Seq_3_cell_metadata.csv"

all_cell_meta <- read.csv(all_cell_path, row.names =1)

fpkm_path <- "output/transcript_fpkm_matrix.csv"
fpkm_mat <- read.table(fpkm_path, header = T)

tpm_path <- "output/transcript_tpm_matrix.csv"
tpm_mat <- read.table(tpm_path, header = T)

sub_organoid <- function(expr_mat, expr_path){
	# browser()
	colnames(expr_mat) <- gsub(".*_S", "S", colnames(expr_mat))	
	
	org_expr_mat <- expr_mat[,!colnames(expr_mat) %in% fetal_cell_meta$Sample_ID]
	fetal_expr_mat <- expr_mat[,c("TRANSCRIPT_ID", fetal_cell_meta$Sample_ID)] 
	
	print(dim(org_expr_mat))
	print(dim(fetal_expr_mat))
	
	org_expr_path <- gsub(".csv", "_organoid.csv", expr_path)
	fetal_expr_path <- gsub(".csv", "_fetal.csv", expr_path)
	
	print(org_expr_path)
	print(fetal_expr_path)
	
	write.table(org_expr_mat, org_expr_path)
	write.table(fetal_expr_mat, fetal_expr_path)
	
}

sub_organoid(fpkm_mat, fpkm_path)
sub_organoid(tpm_mat, tpm_path)

?trace

