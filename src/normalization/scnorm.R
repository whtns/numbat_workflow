#scnorm check read count/depth relationship on HSMM

Conditions = rep(c(1), each= 312)
checkCountDepth(Data = as.matrix(exprs(HSMM)), Condition = Conditions, OutputName = "./single_cell_pipeline/output/scnorm/sunhye_census_norm_check",
                FilterCellProportion = .1, FilterExpression = 2)
stringtie_gene_matrix <-na.omit(stringtie_gene_matrix)
DataNorm <- SCnorm(stringtie_transcripts_matrix, Conditions, OutputName = "./single_cell_pipeline/output/scnorm/sunhye_genes_scnorm",
                   FilterCellNum = 10)
write.table(DataNorm$NormalizedData, file = "./single_cell_pipeline/output/scnorm/sunhye_genes_scnorm.csv")
