#!/usr/bin/Rscript

library(tidyr)
library(dplyr)

#DS
# cell_settings = "~/single_cell_tools/dshayler_input/2_seq_3dformat_050418.csv"
# plot_settings = "~/single_cell_tools/dshayler_input/030618_3d_PCA_No_Bad_Reads_Color_by_age.txt"

#SHL
# cell_settings <- "~/single_cell_tools/SHL_supplied_20180221_FACS_20170407_SHL_input_files/New_cells_sets.csv"
# plot_settings <- "~/single_cell_tools/SHL_supplied_20180221_FACS_20170407_SHL_input_files/New_plot_settings_2d.csv"
# cell_settings <- "~/single_cell_tools/FACS_0407_2017_SHL_input_files/cell_sets_0407_SHL_20180523.csv"
cell_settings <- "~/single_cell_pipeline/tmp/cell_settings_20180618.csv"
plot_settings <- "~/single_cell_pipeline/tmp/plot_settings_20180618.csv"



find_remove_cells <- function(plot_settings, annotation){
  
  test <- readLines(plot_settings)
  
  if (!grepl('remove', test)){
    return(NULL)
  }
  
  vecs <- list()
  mtnames <- c()
  for (i in test){
    if (!grepl("#", i) & grepl('remove', i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[2]]
      vecs <- append(vecs, lbline)
      mtnames <- append(mtnames, d)    
    }
  }
  
  pfx <- tolower(gsub("_.*", "", mtnames))
  valid_pfx  <- which(pfx %in% colnames(annotation))
  
  if (length(valid_pfx) == 0){
    return(NULL)
  }
  
  pfx <- pfx[valid_pfx]
  
  
  
  sfx <- gsub(".*_", "", mtnames[valid_pfx])
  
  remove_cells <- purrr::map2(pfx, sfx , function(x, y) annotation[annotation[tolower(x)] == y,])
  
  
  
  remove_cells <- dplyr::bind_rows(remove_cells)
  
  ind <- apply(remove_cells, 1, function(x) all(is.na(x)))
  remove_cells <- remove_cells[ !ind, ]
  
  remove_cells <- unique(remove_cells[,1])
  
}


match_cols <- function(match_vecs, sv_name){
  out=NULL
  for (i in match_vecs){
    vetor <- i
    vetor <- vetor[vetor != ""]
    key <- data.frame(sample_id=vetor[-1], sv_name=rep(gsub(".*_","", vetor[[1]]), (length(vetor)-1)))  
    out <- rbind(out, key)
  }  
  colnames(out) <- c("sample_id", sv_name) 
  return(out)
}

convert_mt_setting <- function(cell_settings, plot_settings){

  test <- readLines(cell_settings)
  
  vecs <- list()
  mtnames <- c()
  for (i in test){
    if (!grepl("#", i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[1]]
      vecs <- append(vecs, lbline)
      mtnames <- append(mtnames, d)    
    }
    
  }
  
  pfx <- unique(gsub("_.*", "", mtnames[grep("_", mtnames)]))
  pfx <- paste0(pfx, "_")
  test <- list()
  for (i in pfx){
    test1 <- list(which(startsWith(mtnames, i)))
    names(test1) = tolower(gsub("_", "", i))
    test <- append(test, test1)
  }
  
  sub_vecs <- list()
  vec_names <- list()
  for (i in test){
    test_vec <- vecs[i]
    sub_vecs <- append(sub_vecs, list(test_vec))
  }
  
  # names(sub_vecs) <- vec_names[1:length(sub_vecs)]
  names(sub_vecs) <- names(test)
  
  sub_vecs <- sub_vecs[unlist(lapply(sub_vecs, length) != 0)]
  
  param_dfs <- purrr::map2(sub_vecs, names(sub_vecs), match_cols)
  
  
  
  if (is.list(param_dfs) & length(param_dfs) != 0) {
    param_dfs <- param_dfs %>%
      Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="sample_id"), .) %>% 
      arrange(sample_id)  
    
    dup_cells <- which(duplicated(param_dfs[,1]))
    if (any(dup_cells)){
      print(paste0("cells ", paste(param_dfs$sample_id[dup_cells], collapse = " "), " found duplicated in cell sets! They will be removed from analysis"))
      param_dfs <- param_dfs[-dup_cells,]
    }
    
    rownames(param_dfs) <- param_dfs[,1]
    
  }
  
  remove_cells <- find_remove_cells(plot_settings, param_dfs)

  
  return(list("annotation" = param_dfs, "removed_cells" = remove_cells))
  
}

annotation <- convert_mt_setting(cell_settings, plot_settings)


