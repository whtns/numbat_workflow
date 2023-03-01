#!/usr/local/bin/Rscript


#=====================================
library(optparse)

# os exosome --------------------------------------------------------------
# default_cell_info = "~/os_exosome_pipeline/data/Phase1_OS_controls_CHL_HLOH_sent_SBI.csv"
# default_feature = "gene"
# default_groupfile = "~/os_exosome_pipeline/bin/diffex_group_comparison.txt"
# #after filtering
# # default_infile = "~/os_exosome_pipeline/output/gene_count_matrix.csv"
# # before duplicate filtering
# default_infile = "~/os_exosome_pipeline/output/maverix_bams/all_normal_vs_os_bams/gene_count_matrix.csv"
# default_out = "~/os_exosome_pipeline/output/"
# # raw gene counts
# input_rda <- "~/os_exosome_pipeline/output/os_exosome_before_duplicate_filtering_gene_differential_expression_input.rda"
# # raw transcript counts
# # input_rda <- "~/os_exosome_pipeline/output/os_exosome_transcript_differential_expression_input.rda"


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

# DS 20181001------------------------------------------------------------
# default_cell_info = "~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/FACS_20181001_dshayler_Organoid_proj/10_2018_Seq_3_all_cell_metadata.csv"
# default_feature = "transcript"
# default_groupfile = "~/single_cell_projects/quicklinks/FACS_20181001_dshayler_Organoid_proj/output/differential_expression/diffex_groups_20181017.csv"
# default_infile = "~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/FACS_20181001_dshayler_Organoid_proj/output/transcript_count_matrix.csv"
# # default_infile = "~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/FACS_20181001_dshayler_Organoid_proj/output/gene_count_matrix.csv"
# default_out = "~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/FACS_20181001_dshayler_Organoid_proj/output/differential_expression"
# default_cell_settings <- "~/single_cell_tools/dshayler_input/2_seq_3dformat_050418.csv"

# DS 3seq -----------------------------------------------------------------

# default_infile =  "~/single_cell_projects/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20170407_20171031_20181001_dshayler_H_sapiens_merge_proj/output/All_cells_stringtie_transcripts_raw_counts_merged.csv"
# default_cell_info = "~/single_cell_projects/quicklinks/3_seq_dshayler_proj/3_seq_all_cells_metadata.csv"
# default_groupfile = "~/single_cell_tools/dshayler_input/org_FR/diff_expr/101718_differential_expression_rod_cone_FR_org_clusters.txt"
# default_out = "~/single_cell_tools/dshayler_input/org_FR/diff_expr/101818_rod_cone_all_cluster_DE"
# default_cell_settings <- "~/single_cell_tools/dshayler_input/2_seq_3dformat_050418.csv"


# DS 2seq (yellow vs blue groups)-----------------------------------------

default_infile =  "~/single_cell_projects/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20170407_20171031_dshayler_H_sapiens_merge_proj/stringtie_transcripts_raw_counts_2_seq.csv"
default_cell_info = "~/single_cell_projects/quicklinks/2_seq_dshayler_proj/2_seq_meta_111517.csv"
default_groupfile = "~/single_cell_projects/quicklinks/2_seq_dshayler_proj/2_seq_differential_expression_blue_yellow_FR_groupfile.txt"
default_out = "~/single_cell_projects/quicklinks/2_seq_dshayler_proj/output/differential_expression/blue_yellow_FR"
default_cell_settings <- "~/single_cell_projects/quicklinks/2_seq_dshayler_proj/2_seq_differential_exp_yellow_blue_FR_cell_settings.csv"


# DS 2seq (old)------------------------------------------------------------
# default_cell_info = "/home/skevin/tmp/2_seq_meta_111517.csv"
# # default_cell_info = "~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/2_seq_dshayler/2_seq_meta_111517.csv"
# default_feature = "transcript"
# default_groupfile = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_diffex_comparison_file.csv"
# default_infile = "/home/skevin/tmp/FACS_20170407_20171031_dshayler_transcripts_raw_counts.csv"
# # default_infile = "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_transcripts_raw_counts.csv"
# default_out = "~/single_cell_pipeline/output/merged_analyses/"
# # # raw transcript counts
# input_rda <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_20171031_dshayler_differential_expression_input.rda"


# SHL RD 20181009-----------------------------------------------------------------

#default_infile =  "~/single_cell_projects/sc_RB_devel/RD_20180109_SHL_H_sapiens_RB_31_proj/output/transcript_count_matrix.csv"
#default_cell_info = "~/single_cell_projects/sc_RB_devel/RD_20180109_SHL_H_sapiens_RB_31_proj/RD_Seq_sample_info_20180109.csv"
#default_groupfile = "~/single_cell_projects/sc_RB_devel/RD_20180109_SHL_H_sapiens_RB_31_proj/group_settings.tsv"
#default_out = "~/single_cell_projects/sc_RB_devel/RD_20180109_SHL_H_sapiens_RB_31_proj/output/differential_expression/"
#default_feature = "transcript"
#default_cell_settings <- "~/single_cell_projects/sc_RB_devel/RD_20180109_SHL_H_sapiens_RB_31_proj/cell_settings.txt"


# SHL 20171031-----------------------------------------------------------------

# default_infile =  "~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/transcript_count_matrix.csv"
# default_cell_info = "~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/FACS_20171031_sunlee_sample_sheet.csv"
# default_groupfile = "~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/scde/FACS_1031_sunlee_diffex_groups_20181030_LatePT.csv"
# default_out = "~/single_cell_projects/quicklinks/FACS_20171031_sunlee_H_sapiens_proj/output/"
# default_feature = "transcript"
# default_cell_settings <- "cell_settings_2nd_Exp_Late_PT1_for_scde_analysis_20181020.tsv, cell_settings_2nd_Exp_Late_PT2_for_scde_analysis_20181030.tsv"

# # SHL 20170407 SCDE------------------------------------------------------------
# default_infile =  "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/transcript_count_matrix.csv"
# default_cell_info = "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/FACS_20170407_sunlee_sample_sheet.csv"
# default_groupfile = "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/scde_lateshctrl_shRb1_scde_groupfile.csv"
# default_out = "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/scde_late_Ctrl"
# default_cell_settings <- "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/cell_settings_1st_exp_latectrl_PT1_scde.csv"


# SHL 20170407 scde-----------------------------------------------------------------

# default_infile =  "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/transcript_count_matrix.csv"
# default_cell_info = "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/FACS_20170407_sunlee_sample_sheet.csv"
# default_groupfile = "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/FACS_0407_sunlee_diffex_groups_20181025_latePT_edgeR.csv"
# default_out = "~/single_cell_projects/quicklinks/FACS_20170407_sunlee_H_sapiens_proj/output/"
# default_feature = "transcript"
# default_cell_settings <- "cell_settings_1st_Exp_PT1_for_scde_analysis_20181002.tsv,cell_settings_1st_Exp_PT2_for_scde_analysis_20181002.tsv"


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
  make_option(c("-s", "--cell_settings"), type="character", default=default_cell_settings,
  						help="tab delimited cell settings file [default= %default]", metavar="character"),
  make_option(c("-f", "--feature"), type="character", default="transcript", 
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
suppressMessages(library(cataract))
# suppressMessages(library(sva))
edb <- EnsDb.Hsapiens.v86
suppressMessages(library(zinbwave))


# load functions -------------------------------------------------

plot_heatmap <- function(filt_vsd, heatmap_path, file_tag){
  heatmap_path <- paste0(heatmap_path, file_tag, "_heatmap.html")
  heatmaply(filt_vsd, scale = "row", seriate = "none", dendrogram = "both", 
            hclust_method = "ward", file = heatmap_path, 
            colors = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255) )
  
}


find_top_var_genes <- function(vsd_obj, top_n){
  vsd <- head( order( rowVars( assay(vsd_obj) ), decreasing=TRUE ), top_n )
  vsd <- assay(vsd_obj)[vsd,]
  syms <- lookup_genes(rownames(vsd))
  rownames(vsd) <- paste0(syms, "\n", rownames(vsd))
  return(vsd)
} 

find_top_sig_genes <- function(vsd_obj, res, top_n){
  vsp <- rownames(head(res[order(res$padj),], top_n))
  vsp <- assay(vsd_obj)[vsp,]
  syms <- lookup_genes(rownames(vsp))
  rownames(vsp) <- paste0(syms, "\n", rownames(vsp))
  return(vsp)
} 

find_user_supplied_features <- function(vsd_obj, in_features, top_n){
  in_features <- readr::read_csv(in_features, header = FALSE)
  vsp <- in_features[[1]]
  vsp <- assay(vsd_obj)[vsp,]
  syms <- lookup_genes(rownames(vsp))
  rownames(vsp) <- paste0(syms, "\n", rownames(vsp))
  return(vsp)
} 

take_input <- function(prompt = question, interactive = F){
  if(interactive){
    param <- readline(prompt=question)
  } else {
    param <- readLines("stdin", n=1)
  }
  return(param)
}

# load data ---------------------------------------------------
suppressWarnings(dir.create(opt$out))

# read in group_settings file
diffex_settings <- readr::read_tsv(opt$groupfile)
# convert dataframe into lists

# convert group_settings dataframe into lists
diffex_settings <- purrr::transpose(diffex_settings)

cell_settings <- unlist(strsplit(opt$cell_settings, ","))

# append cell_settings info onto diffex_settings
diffex_settings <- purrr::map2(diffex_settings, cell_settings, append_mt_settings)

# read in expression/raw count matrix
counts <- readr::read_csv(opt$infile)

# read in cell metadata
sample_sheet <- readr::read_csv(opt$cell_info)

# trace("prep_counts_and_sample_sheet", quote(browser(skipCalls = 4)),
#  	exit = quote(browser(skipCalls = 4)))

diffex_in <- purrr::map(diffex_settings, prep_counts_and_sample_sheet, counts, sample_sheet)


#specify model design
#=============================

# trace("run_edger_qlf", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))
# 
# trace("lookup_symbols_from_genes", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))



while (TRUE) {
  question = "Choose from following:\n
            [DLRT] Run DESEQ2 Likelihood Ratio Test\n
            [SCDE] Run SCDE pairwise comparison\n
            [DC] Run DESEQ2 Pair Comparison\n
            [EC] Run EdgeR Pair Comparison\n
						[ECZ] Run EdgeR Pair comparison with ZINB correction\n
            [EQLF] Run EdgeR Quasi-Likelihood F Test\n
            [LTF] Run limma trend F test\n
            [LVF] Run limma voom F test\n
            [HM] Create Heatmaps\n
            [X] Exit\n "
  # if (ctr < 3){
  #   print(question)
  # }
  
  cat(question)
  # non-interactive (in terminal)
  action <- toupper(take_input(prompt=question))
  # interactive (in rstudio)
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
    
    
  } else if (action == "DC"){
      
      check_file_exist <- function(destfile){
  
        if(file.exists(destfile)){
          dds <- readRDS(destfile)
        } else {
          cat("preparing deseq2 input\n")
      	}
      }
      
      outdir <- paste0(opt$out, "/DESEQ2/", opt$feature, "/")
      dir.create(outdir,showWarnings = FALSE)
      
      ddsRDS <- paste0(outdir, user_groups$analysis_id, "_dds.rds")
  
      ddsl <- lapply(ddsRDS, check_file_exist)
  

      if(any(sapply(ddsl, is.null))){
        cat("running deseq2 for between group comparison")
        ddsl <- purrr::map2(user_groups$full_model, summarized_exp, function(x,y){dds <- DESeqDataSetFromMatrix(countData = y$counts,
                                                                                                                colData = y$sample_sheet,
                                                                                                                design = as.formula(x)); return(dds)})
        
        names(ddsl) <- user_groups$analysis_id
        
        # minimal filtering to remove outlier genes present in few cells ----------
        ddsl <- lapply(ddsl, function(x) x[cataract::meetsCondition(counts(x), gtet=10, nCells=10),])
        
        run_DESeq <- function(dds_obj){
          dds <- DESeq(dds_obj, parallel = TRUE, BPPARAM = MulticoreParam(), betaPrior = FALSE)
          return(dds)
        }
        
        ddsl <- purrr::map(ddsl, run_DESeq)
        outdir <- paste0(opt$out, "/DESEQ2/", opt$feature, "/")
        dir.create(outdir,showWarnings = FALSE)
        ddsRDS <- paste0(outdir, user_groups$analysis_id, "_dds.rds")
        purrr::map2(ddsl, ddsRDS, saveRDS)
        ddsl <- lapply(ddsRDS, readRDS)
        
      }
  
      cell_info <- cataract::safe_read(opt$cell_info, row.names = 1)
  
      results_models <- gsub("~", "", user_groups$full_model)
      results_models <- mapply(c, "full_model" = results_models, "analysis_id" = user_groups$analysis_id, SIMPLIFY = FALSE)
      
      diffex_comps <- purrr::map2(ddsl, results_models, cataract::make_all_deseq2_comparisons, cell_info, opt$feature) 
            
  }  else if (action == "ECZ"){

    question = "prefiltering minimum count per gene (leave blank for default: 1) "
    cat(question)
    min.count <- take_input(prompt=question)
    if(min.count == ""){
      min.count = 1
    } else {
      min.count = as.numeric(min.count)
    }
    
    question = "prefiltering number of cells with minimum count (leave blank for default: 2) "
    cat(question)
    min.total.count <- take_input(prompt=question)
    if(min.total.count == ""){
      min.total.count = 2
    } else {
      min.total.count = as.numeric(min.total.count)
    }
    
    tic("finished EdgeR pairwise comparison test with")
    # cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
    
    qlfl <- purrr::map(diffex_in, run_edger_pair_zinb, 
    									 min.count, min.total.count, feature = opt$feature, outdir = opt$out)
    
    toc()
      
      
  } else if (action == "EC"){
  	
  	question = "prefiltering minimum count per gene (leave blank for default: 1) "
  	cat(question)
  	min.count <- take_input(prompt=question)
  	if(min.count == ""){
  		min.count = 1
  	} else {
  		min.count = as.numeric(min.count)
  	}
  	
  	question = "prefiltering number of cells with minimum count (leave blank for default: 2) "
  	cat(question)
  	min.total.count <- take_input(prompt=question)
  	if(min.total.count == ""){
  		min.total.count = 2
  	} else {
  		min.total.count = as.numeric(min.total.count)
  	}
  	
  	tic("finished EdgeR pairwise comparison test with")
  	# cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
  	
  	# model_comparison <- purrr::map2(user_groups$full_model, user_groups$pair_comp, c)
  	
  	# qlfl <- purrr::map2(summarized_exp, model_comparison, run_edger_pair_comparison, 
  	#                     min.count, min.total.count, feature = opt$feature)
  	
  	
  	qlfl <- purrr::map(diffex_in, run_edger_pair_comparison, 
  											min.count, min.total.count, feature = opt$feature, outdir = opt$out)
  	
  	toc()
  	
  
  	
  } else if (action == "SCDE"){
  	
  	question = "prefiltering minimum count per gene (leave blank for default: 1) "
  	cat(question)
  	min.count <- take_input(prompt=question)
  	if(min.count == ""){
  		min.count = 1
  	} else {
  		min.count = as.numeric(min.count)
  	}
  	
  	question = "prefiltering number of cells with minimum count (leave blank for default: 2) "
  	cat(question)
  	min.total.count <- take_input(prompt=question)
  	if(min.total.count == ""){
  		min.total.count = 2
  	} else {
  		min.total.count = as.numeric(min.total.count)
  	}
  	
  	tic("finished SCDE pairwise comparison test with")
  	# cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
  	
  	# model_comparison <- purrr::map2(user_groups$full_model, user_groups$pair_comp, c)
  	
  	# qlfl <- purrr::map2(summarized_exp, model_comparison, run_edger_pair_comparison, 
  	#                     min.count, min.total.count, feature = opt$feature)
  	
  	
  	ediff <- purrr::map(diffex_in, run_scde, 
  										 min.count, min.total.count, feature = opt$feature, outdir = opt$out)
  	
  	toc()
  	
  	
  } else if (action == "EQLF"){
      question = "prefiltering minimum count per gene (leave blank for default: 1) "
      cat(question)
      min.count <- take_input(prompt=question)
      if(min.count == ""){
        min.count = 1
      } else {
        min.count = as.numeric(min.count)
      }
      
      question = "prefiltering number of cells with minimum count (leave blank for default: 2) "
      cat(question)
      min.total.count <- take_input(prompt=question)
      if(min.total.count == ""){
        min.total.count = 2
      } else {
        min.total.count = as.numeric(min.total.count)
      }
      
      tic("finished EdgeR LRT test with")
      # cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
      
      
      outdir <- paste0(opt$out, "/EDGER/", opt$feature, "/")
      dir.create(outdir,showWarnings = FALSE)
      
      
      print(paste("saving output csvs to", qlfcsv))
      
      purrr::map2(qlfl, qlfcsv, function(x,y) write.csv(x, y))
      toc()
    
  } else if (action == "LTF"){
    
    question = "prefiltering minimum count per gene for at least some samples (leave blank for default: 10) "
    cat(question)
    min.count <- take_input(prompt=question)
    if(min.count == ""){
      min.count = 1
    } else {
      min.count = as.numeric(min.count)
    }
    
    question = "prefiltering minimum count per gene for sum of all samples (leave blank for default: 15) "
    cat(question)
    min.total.count <- take_input(prompt=question)
    if(min.total.count == ""){
      min.total.count = 2
    } else {
      min.total.count = as.numeric(min.total.count)
    }
    
    tic("finished limma F-test with")
    # cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
    
    
    
    outdir <- paste0(opt$out, "/LIMMA/", opt$feature, "/")
    dir.create(outdir,showWarnings = FALSE)
    
    
    
    purrr::map2(tops, topcsv, function(x,y) write.csv(x, y))
    toc()
    
  } else if (action == "LVF"){
    
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
    
    tic("finished limma F-test with")
    # cat(paste0("running EdgeR LRT test using full model: ", full_model_params, " and reduced model: ", reduced_model_params, "\n"))
    
    
    doutdir <- paste0(opt$out, "/LIMMA/", opt$feature, "/")
    dir.create(outdir,showWarnings = FALSE)
    
    topcsv <- paste0(outdir, user_groups$analysis_id, "_limma_voom_f_test.csv")
    print(paste("saving output csvs to", topcsv))
    
    purrr::map2(tops, topcsv, function(x,y) write.csv(x, y))
    toc()
    
  } else if (action == "HM"){
    
    question = "Select method of filtering genes: variance (var), signficance (sig), user-supplied (user) "
    cat(question)
    filter_method <- take_input(prompt=question)
    if(filter_method == "user"){
      filter_method == "var"
    } 
    
    question = "select differential expression method: deseq2, edgeR, limma-trend, limma-voom "
    cat(question)
    diffex_method <- take_input(prompt=question)

    vsdl <- lapply(ddsl, vst, blind = F)

    heatmap_paths <- paste0(opt$out, user_groups$analysis_id)
        
    if (filter_method == "var"){
      topVarGenes <- lapply(vsdl, find_top_var_genes, 500)
      names(topVarGenes) <- gsub("_res.rds", "_top_var", basename(resRDS))
      purrr::map2(topVarGenes, heatmap_paths, plot_heatmap, "_top_var_genes")
    } else if (filter_method == "sig"){
      topSigGenes <- purrr::map2(vsdl, resl, find_top_sig_genes, 500)
      names(topVarGenes) <- gsub("_res.rds", "_top_sig", basename(resRDS))
      purrr::map2(topSigGenes, heatmap_paths, plot_heatmap, "_top_sig_genes")
    } else if (filter_method == "user"){
      user_genes <- "~/single_cell_tools/test_pcs_4_6_7.csv"
      topUserGenes <- purrr::map(vsdl, find_user_supplied_features, user_genes, 500)
      names(topUserGenes) <- gsub("_res.rds", "_top_user", basename(resRDS)) 
      purrr::map2(topUserGenes, heatmap_paths, plot_heatmap, "_top_user_genes")
    }
  } else if (action == "I"){
  browser()
  print(a)
  }
}

# trace("run_limma_voom", quote(browser(skipCalls = 4)),
# 			exit = quote(browser(skipCalls = 4)))
# 
# trace("make_all_deseq2_comparisons", quote(browser(skipCalls = 4)),
# exit = quote(browser(skipCalls = 4)))
# 
# trace("run_DESeq_LRT", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))
# # #
# trace("output_deseq2_results", quote(browser(skipCalls = 4)),
# exit = quote(browser(skipCalls = 4)))
# 
