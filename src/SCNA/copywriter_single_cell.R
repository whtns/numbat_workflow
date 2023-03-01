#!/usr/bin/Rscript

library(stchlk.copywriter)
library(optparse)


default_binsize = 1010000
default_inputdir = "~/single_cell_pipeline/output/GT_20180109_SHL_H_sapiens_RB_31_output/"
default_bamfileregex = "*under_3.bam$"


#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-b", "--binsize"), type="character", default=default_binsize,
              help="bin size [default= %default]", metavar="character"),
  make_option(c("-i", "--inputdir"), type="character", default=default_inputdir,
              help="input directory for bam to be analyzed [default= %default]", metavar="character"),
  make_option(c("-r", "--bamfileregex"), type="character", default=default_bamfileregex,
              help="regex that defines bams to be analyzed [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

outputdir = paste0(opt$inputdir, "/copywriter/", humanreadable(as.numeric(opt$binsize)))

# trace("run_copywriter", quote(browser(skipCalls = 4)),
#       exit = quote(browser(skipCalls = 4)))

run_copywriter(opt$binsize, opt$inputdir, outputdir, opt$bamfileregex, ref_genome = "hg38")
# 
# untrace("run_copywriter")
