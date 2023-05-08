source("packages.R")
source("functions.R")
library(targets)
tar_load("large_filter_expressions")
tar_load("cluster_dictionary")
# debug(find_diffex_clones)
# debug(make_clone_comparison)

debug(make_numbat_plot_files)

make_numbat_plot_files("output/numbat_sridhar/SRR14800541/done.txt", cluster_dictionary, large_filter_expressions, extension = "_filtered")
