#!/usr/bin/Rscript

library(optparse)

# default Dominic ---------------------------------------------------------
default_infile = "~/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/2_seq_dshayler/stringtie_transcripts_raw_counts_2_seq.csv"
default_cell_info = "~/single_cell_pipeline/output/merged_analyses/peleg_dominic_df_analysis4.csv"
default_groupfile = "~/single_cell_pipeline/scde_input/scde_peleg_analysis1_20180228.csv"
default_out = "/home/pwiner/"
default_feature = "transcript"

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
  make_option(c("-f", "--feature"), type="character", default=default_feature, 
              help="the feature (gene or transcript) under analysis", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

library(cataract)

print(packageVersion("flexmix"))

user_groups <- read.table(textConnection(gsub(" |,", "\t", readLines(opt$groupfile))), sep="\t", header = TRUE, stringsAsFactors = FALSE, fill=TRUE)
comparisons <- user_groups$comparison
analysis_ids = user_groups$analysis_id
rds_files <- list.files(path = paste0(opt$out, "/diffex_by_trs"), pattern=".*stringtie.rds$", full.names = TRUE)

rds_files <- sapply(analysis_ids, function(x) rds_files[grep(x, rds_files)])

rds_names <- basename(rds_files)
rds_names <- gsub("_stringtie.rds", "", rds_names)

scde_preps <- purrr::map(rds_files, readRDS)
scde_preps <- purrr::map2(scde_preps, rds_names, function(x,y){c(x, analysis_id=y)})

cd <- scde_preps$analysis1$cd

# knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 1, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)

saveRDS(knn, "knn.rds")

knn <- readRDS("knn.rds")

varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = TRUE)
top_overd_trx <- sort(varinfo$arv, decreasing = TRUE)[1:10]

top_overd_df <- data.frame("genes" = lookup_genes(names(top_overd_trx)), "overdispersion" = top_overd_trx)


# Evaluate overdispersion of 'de novo' gene sets --------------------------

clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 1, plot = TRUE)