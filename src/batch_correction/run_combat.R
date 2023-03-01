#!/usr/local/bin/Rscript

#load required libraries and functions
#=====================================
library(optparse)

# os exosome --------------------------------------------------------------
# default_cell_info = "~/os_exosome_pipeline/data/Phase1_OS_controls_CHL_HLOH_sent_SBI.csv"
# default_feature = "gene"
# default_groupfile = "~/os_exosome_pipeline/bin/diffex_group_comparison.txt"
# default_infile = "~/os_exosome_pipeline/output/gene_count_matrix.csv"
# default_out = "~/os_exosome_pipeline/output/"
# raw gene counts
# input_rda <- "~/os_exosome_pipeline/output/os_exosome_gene_differential_expression_input.rda"
# raw transcript counts
# input_rda <- "~/os_exosome_pipeline/output/os_exosome_transcript_differential_expression_input.rda"


# SHL 20170407 ------------------------------------------------------------
# default_cell_info = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/cell_info.csv"
# default_feature = "transcript"
# default_groupfile = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/DESEQ2/shl_deseq2_comparison_file.csv"
# default_infile = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/transcript_count_matrix.csv"
# # default_infile = "/home/skevin/tmp/FACS_20170407_sunlee_census_matrix_mod.csv" # census
# # default_infile = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/gene_count_matrix.csv"
# default_out = "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/"
# # raw transcript counts
# input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_differential_expression_input.rda"
# # census transcript counts
# # input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_census_differential_expression_input.rda"
# # raw gene counts
# # input_rda <- "~/single_cell_pipeline/output/FACS_20170407_sunlee_H_sapiens_output/shl_0407_genes_differential_expression_input.rda"

# DS 2seq ------------------------------------------------------------
default_cell_info = "/home/skevin/tmp/2_seq_meta_111517.csv"
# default_cell_info = "~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/2_seq_dshayler/2_seq_meta_111517.csv"
default_feature = "transcript"
default_groupfile = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_diffex_comparison_file.csv"
default_infile = "/home/skevin/tmp/FACS_20170407_20171031_dshayler_transcripts_raw_counts.csv"
# default_infile = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_transcripts_raw_counts.csv"
default_out = "~/single_cell_pipeline/output/merged_analyses/"
# # raw transcript counts
input_rda <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_differential_expression_input.rda"


# diffex_input <- mget(ls(pattern = "default"))
# save(diffex_input, file = input_rda)

if (file.exists(input_rda)){
	load( input_rda)
	list2env(diffex_input, globalenv())
} else {
	default_expr_mat = default_annotation =  default_cell_info =  default_plot_settings = default_out = NA
	
}


# default input files -----------------------------------------------------

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
	make_option(c("-i", "--infile"), type="character", default=default_infile,
							help="gene expression input filename [default= %default]", metavar="character"),
	make_option(c("-c", "--cell_info"), type="character", default=default_cell_info,
							help="metadata about cells in input file [default= %default]", metavar="character"),
	make_option(c("-g", "--groupfile"), type="character", default=default_groupfile,
							help="file that defines comparison groups [default= %default]", metavar="character"),
	make_option(c("-o", "--out"), type="character", default=default_out,
							help="output file name [default= %default]", metavar="character"),
	make_option(c("-f", "--feature"), type="character", default=NA, 
							help="the feature (gene or transcript) under analysis", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
	print_help(opt_parser)
	missing_args <- paste(names(which(sapply(opt, is.na))), collapse = " ")
	stop(paste0("Please provide all necessary arguments. Missing: ", missing_args), call.=FALSE)
}


# load required libraries -------------------------------------------------
print("loading required libraries")
suppressMessages(library(BiocParallel))
suppressMessages(library(tidyverse))
suppressMessages(library(gtools))
suppressMessages(library(cataract))
suppressMessages(library(data.table))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(gplots))
suppressMessages(library(tictoc))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(Matrix))
suppressMessages(library(sva))
edb <- EnsDb.Hsapiens.v86
suppressMessages(library(zinbwave))


# load required functions -------------------------------------------------
plot_pca <- function(counts_df, sample_sheet, color_var = color_var, pca_title = "PCA"){
  # browser()
  PC <- prcomp(t(counts_df))$x
  
  sample_seq_df <- sample_sheet[c("sample_id", color_var)]
  sample_seq_df[[color_var]]<-as.factor(sample_seq_df[[color_var]])
  
  rownames(sample_seq_df) <- sample_seq_df[,"sample_id"]
  
  sample_seq_df <- sample_seq_df[color_var]
  
  PC_seq <- merge(PC, sample_seq_df, by = 0)
  
  g <- ggplot(PC_seq,aes_string(x="PC1",y="PC2",col=color_var))+
    geom_point(size=3,alpha=0.5)+ #Size and alpha just for fun
    # scale_color_manual(values = c("#FF1BB3","#A7FF5B"))+ #your colors here
    theme_classic() + 
    ggtitle(pca_title) + 
    NULL
  print(g)
}

run_combat <- function(summarized_exp, full_model, min.count, min.total.count, color_var = "Seq_Number", feature = "transcript"){
	# browser()
	full_model <- gsub("~", "", full_model)
	dgelist <- DGEList(counts=summarized_exp$counts, group=summarized_exp$sample_sheet[[full_model]])
	
	design <- model.matrix(~group, data=dgelist$samples)
	colnames(design) <- levels(dgelist$samples$group)
	
	keep <- filterByExpr(dgelist, design, min.count = min.count, min.total.count = min.total.count)
	dgelist <- dgelist[keep,keep.lib.sizes=FALSE]
	
	dgelist <- calcNormFactors(dgelist)
	
	v <- voom(dgelist, design, plot=TRUE)
	
	batch = summarized_exp$sample_sheet$Seq_Number
	
	modcombat = model.matrix(~Fetal_Age, data=summarized_exp$sample_sheet)
	
	combat_edata = ComBat(dat=as.matrix(v$E), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

	raw_pca_seq <- plot_pca(v$E, summarized_exp$sample_sheet, color_var="Seq_Number", pca_title = "raw data")
	combat_pca_seq<-plot_pca(combat_edata, summarized_exp$sample_sheet, color_var = "Seq_Number", pca_title = "combat adjusted")
	
	raw_pca_age <- plot_pca(v$E, summarized_exp$sample_sheet, color_var="Fetal_Age", pca_title = "raw data")
	combat_pca_age<-plot_pca(combat_edata, summarized_exp$sample_sheet, color_var = "Fetal_Age", pca_title = "combat adjusted")
	
		
	return(list("raw" = v$E, "combat" = combat_edata))	
}

lookup_symbols_from_genes <- function(geneid){
	
	geneid_filt <- str_replace(geneid, "\\.\\d+.*", "")
	syms <- genes(edb, filter = GeneIdFilter(geneid_filt), columns = c("symbol"))
	syms <- as.data.frame(syms)
	syms <- data.frame(row.names = syms[["gene_id"]], "symbol" = syms[["symbol"]])
	
	return(syms)
}

lookup_symbols_from_transcripts <- function(txid){
	
	txid <- str_replace(txid, "\\.\\d+.*", "")
	syms <- transcripts(edb, filter = TxIdFilter(txid), columns = c("symbol"))
	syms <- as.data.frame(syms)
	syms <- data.frame(row.names = syms[["tx_id"]], "symbol" = syms[["symbol"]])
	
	return(syms)
}

take_input <- function(prompt = question, interactive = TRUE){
	if(interactive){
		param <- readline(prompt=question)
	} else {
		param <- readLines("stdin", n=1)
	}
	return(param)
}

#=================================================

# make output direcotry ---------------------------------------------------
suppressWarnings(dir.create(opt$out))

user_groups <- read.table(textConnection(gsub(" ", "\t", readLines(opt$groupfile))), sep="\t", header = TRUE, stringsAsFactors = FALSE, fill=TRUE)

counts <- cataract::safe_read(opt$infile)
# counts <- read.table(opt$infile)

sample_sheet <- read.csv(opt$cell_info, row.names=1) 
# 
# trace("prep_counts_and_sample_sheet", quote(browser(skipCalls = 4)),
# exit = quote(browser(skipCalls = 4)))

summarized_exp <- purrr::map2(user_groups$subset_parameter, user_groups$subset_param_values, prep_counts_and_sample_sheet, counts, sample_sheet, strip_col = TRUE)
# 
# test1 <- summarized_exp[[1]]$counts
# test2 <- summarized_exp[[1]]$sample_sheet
# 
# rse <- SummarizedExperiment(assays=SimpleList(counts=test1), colData=test2)

filter_summarized_exp <- function(sexp, keep_cells){
	browser()
	sexp[[1]]$counts <- sexp[[1]]$counts[,colnames(sexp[[1]]$counts) %in% keep_cells]
	sexp[[1]]$sample_sheet <- sexp[[1]]$sample_sheet[rownames(sexp[[1]]$sample_sheet) %in% keep_cells,]
	return(sexp)
}

summarized_exp_filtered <- filter_summarized_exp(summarized_exp, colnames(filtered_table))

# continue with deseq2 analysis -------------------------------------------

# full_list <- list()
# for (row in 1:nrow(user_groups)){
#     group_list <- as.list(user_groups[row,])
#     names(group_list) <- names(user_groups)
#     print(row)
#     full_list[[row]] <- group_list
#   }  

full_model_params <- gsub("^[[:punct:]]", "", user_groups$full_model)
reduced_model_params <- gsub("^[[:punct:]]", "", user_groups$reduced_model) 


# factorize all model parameters ------------------------------------------
summarized_exp <- purrr::map2(summarized_exp, full_model_params, function(x,y){ x[["sample_sheet"]][[y]] <- as.factor(x[["sample_sheet"]][[y]]); return(x)})

#specify model design
#=============================

# trace("run_edger_qlf", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))
# 
# trace("lookup_symbols_from_genes", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))



while (TRUE) {
	question = "Choose from following:\n
	[C] Run combat\n
	[X] Exit\n "
	# if (ctr < 3){
	#   print(question)
	# }
	
	cat(question)
	action <- toupper(take_input(prompt=question))
	# action <- readline(prompt=question)
	
	if(action == "X"){
		break
	} else if (action == "DLRT"){
		tic("finished LRT test with")
		# cat(paste0("running deseq2 with LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
		# 
		# ddsl <- purrr::map2(summarized_exp, function(x,y){dds <- DESeqDataSetFromMatrix(countData = y$counts,
		#                                                                                                         colData = y$sample_sheet,
		#                                                                                                         design = as.formula(x)); return(dds)})
		# 
		# 
		# names(ddsl) <- user_groups$analysis_id
		# 
		#     # minimal filtering to remove outlier genes present in few cells ----------
		# ddsl <- lapply(ddsl, function(x) x[cataract::meetsCondition(counts(x), gtet=10, nCells=10),])
		
		deseq_args <- lapply(split(user_groups, row.names(user_groups)), unlist)
		
		deseq_args <- purrr::map2(summarized_exp, deseq_args, c)
		
		
		ddsl <- purrr::map(deseq_args, map_DESeq, feature = opt$feature)
		
		toc()
		
		
	} else if (action == "C"){
		
		min.count.default = 30
		
		question = paste0("prefiltering minimum count per gene for at least some samples (leave blank for default: ", min.count.default, ")")
		cat(question)
		min.count <- take_input(prompt=question)
		if(min.count == ""){
			min.count = min.count.default
		} else {
			min.count = as.numeric(min.count)
		}
		
		min.total.count.default = 45
		
		question = paste0("prefiltering minimum count per gene for sum of all samples (leave blank for default: ", min.total.count.default, ")")
		cat(question)
		min.total.count <- take_input(prompt=question)
		if(min.total.count == ""){
			min.total.count = min.total.count.default
		} else {
			min.total.count = as.numeric(min.total.count)
		}
		
		tic("finished combat F-test with")
		# cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
		
		tops <- purrr::map2(summarized_exp, user_groups$full_model, run_combat, min.count, min.total.count, feature = opt$feature)
		
		doutdir <- paste0(opt$out, "/COMBAT/", opt$feature, "/")
		dir.create(outdir,showWarnings = FALSE)
		
		topcsv <- paste0(outdir, user_groups$analysis_id, "_combat.csv")
		print(paste("saving output csvs to", topcsv))
		
		purrr::map2(tops, topcsv, function(x,y) write.csv(x, y))
		toc()
		
	} else if (action == "I"){
		browser()
		print(a)
	}
}
#####-non_filtered combat run
combat_edata <- run_combat(summarized_exp[[1]], user_groups$full_model, min.count=30, min.total.count=45, color_var = "Seq_Number", feature = "gene")
#####-filtered combat run
combat_edata <- run_combat(filtered_table, user_groups$full_model, min.count=30, min.total.count=45, color_var = "Seq_Number", feature = "gene")



purrr::map2(combat_edata, c("raw data", "combat_adjusted"), function(x,y) plot_pca(x, sample_sheet, "Fetal_Age", pca_title = y))


census_path <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_census_matrix.csv"
census_matrix <- read.csv(census_path, row.names = 1)

plot_pca(census_matrix, sample_sheet, "Seq_Number")

raw_pca <- prcomp(t(combat_edata$raw))$x
