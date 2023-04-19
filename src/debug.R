source("packages.R")
source("functions.R")
library(targets)
tar_load("large_in_segment_diffex_clones_for_each_cluster")
tar_load("large_in_segment_diffex_clones")
# debug(find_diffex_clones)
# debug(make_clone_comparison)


make_volcano_diffex_clones(large_in_segment_diffex_clones_for_each_cluster,
                                      "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
                                      large_in_segment_diffex_clones,
                                      "results/diffex_bw_clones_large_in_segment.pdf")
