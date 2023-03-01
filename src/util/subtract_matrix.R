#!/usr/bin/Rscript

#load optparse
#=====================================
library(optparse)

default_infile <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_sunlee_FACS_20171031/FACS_20170407_sunlee_FACS_20171031_sunlee_census_matrix.csv"

default_metafile <- "~/single_cell_pipeline/output/merged_analyses/FACS_20170407_sunlee_FACS_20171031/FACS_20170407_sunlee_FACS_20171031_sunlee_census_matrix_meta.csv"

default_out = "test"
# default input files -----------------------------------------------------

#'  section for parsing command line options when calling script
#'  ###################################
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=default_infile,
              help="census matrix file [default= %default]", metavar="character"),
  make_option(c("-m", "--meta"), type="character", default=default_metafile,
              help="metadata about cells in input file [default= %default]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=default_out,
              help="output directory [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE);

if (any(sapply(opt, is.na))){
  print_help(opt_parser)
  stop("Please provide all necessary arguments.", call.=FALSE)
}

# load required libraries and functions -----------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(gtools))
suppressMessages(library(tidyr))

# load experiment data
print("reading census matrix")
census_matrix <- read.csv(opt$infile)

census_meta <- read.csv(opt$meta)

# correct misassigned rownames
rownames(census_matrix) <- census_matrix[,1]
census_matrix <- census_matrix[,-1]

# convert to matrix type
census_matrix <- as.matrix(census_matrix)

# log transform the census_matrix
census_matrix <- log2(census_matrix + 1)

census_matrix <- data.frame(census_matrix)

# generate dummy data
n <- 1000; p <- 10 # 1000 observations, 20 variables

# generate all data from N(0,1)
set.seed(5)
dat <- as.data.frame(matrix(rnorm(n*p),ncol = p))
ctrl_metadat <- data.frame(Var1 = c("day_0"), Var2 = c("shCtrl"))

treatment_groups <- c("shCtrl", "RBKD")
days <- c("day_4", "day_8")
metadat <- data.frame(expand.grid(days, treatment_groups))

metadat <- rbind(ctrl_metadat, metadat)

split_cols <- split(colnames(dat), ceiling(seq_along(colnames(dat))/(p/nrow(metadat))))

metadat = data.frame(metadat, I(split_cols)) %>% 
  unnest()
colnames(metadat) <- c("day", "treatment_group", "sample_id")

# load required functions
find_mean_diffs <- function(metadat, dat){
  unique_days <- unique(metadat$day)
  unique_days <- unique_days[mixedorder(gsub(".*_", "", unique_days))]
  l <- vector("list", length(unique_days)-1)
  for (i in seq_along(head(unique_days, -1))){
    print(unique_days[i], max.levels = 0)
    ctrl_means_1 <- rowMeans(dat[,metadat[(metadat[["day"]] == unique_days[i] & metadat[["treatment_group"]] == "shCtrl"),]$sample_id])
    ctrl_means_2 <- (rowMeans(dat[,metadat[(metadat[["day"]] == unique_days[i+1] & metadat[["treatment_group"]] == "shCtrl"),]$sample_id]) - ctrl_means_1)
    l[[i]] <- ctrl_means_2
  }
  return(l)
}

subtract_ctrl_diff <- function(metaflag, ctrl_mean_diff, dat, metadat){
  subset_dat <- dat[,metadat[["day"]] == metaflag & metadat[["treatment_group"]] != "shCtrl"]
  corr_subset <- sweep(subset_dat, 1, ctrl_mean_diff, "-")
  return(corr_subset)
  
}

run_ctrl_correction <- function(dat, metadat){

  kd_days <- unique(metadat[metadat[["treatment_group"]] != "shCtrl",]$day)
  ctrl_mean_diffs <- find_mean_diffs(metadat, dat)
  corrected_KD <- purrr::map2_dfc(kd_days, ctrl_mean_diffs, subtract_ctrl_diff, dat, metadat)
  
  rownames(corrected_KD) <- rownames(dat)
  
  positive_log_trans <- corrected_KD + abs(min(corrected_KD))+1
  
  return(positive_log_trans)
}


# run dummy data


# run experiment data
print("subtracting ctrl values for each day")
corrected_census <- run_ctrl_correction(census_matrix, census_meta)


outdir <- opt$outdir
dir.create(outdir)
outfile <- paste0(outdir, "/", gsub(".csv", "_subtracted.csv", basename(opt$infile)))
print(paste0("writing subtracted matrix to: ", outfile))

write.csv(corrected_census, outfile)