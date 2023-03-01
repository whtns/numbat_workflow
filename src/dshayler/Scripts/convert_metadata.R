#!/usr/bin/Rscript

cell_settings <- "~/single_cell_tools/SHL_supplied_20180221_FACS_20170407_SHL_input_files/New_cells_sets.csv"

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

run_conversion <- function(cell_settings){
  test <- readLines(cell_settings)
  
  vecs <- list()
  colnames <- c()
  for (i in test){
    if (!grepl("#", i)){
      lbline = strsplit(i, "\t")
      d = unlist(lbline)[[1]]
      vecs <- append(vecs, lbline)
      colnames <- append(colnames, d)    
    }
    
  }
  
  branch_idx <-  which(grepl("Branch_", colnames))
  day_idx <-  which(grepl("day_", colnames))
  cluster_idx <-  which(grepl("cluster_", colnames))
  treat_idx <- colnames %in% c("sh733", "sh737", "sh842", "shCtrl")
  
  
  branch_vecs <- vecs[branch_idx]
  day_vecs <- vecs[day_idx]
  cluster_vecs <- vecs[cluster_idx]
  treat_vecs <- vecs[treat_idx]
  
  sv_names <- c("branch", "day", "cluster", "treatment_group")
  sub_vecs <- list(branch_vecs, day_vecs, cluster_vecs, treat_vecs)
  names(sub_vecs) <- sv_names
  # vetor <- c(1,2,3)
  # key <- data.frame(vetor=vetor, mat=c('a', 'b', 'c'))
  # data <- data.frame(id=c('a', 'b', 'a', 'c', 'a'))
  
  test <- purrr::map2(sub_vecs, sv_names, match_cols)
  
  test2 <- test %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="sample_id"), .) %>% 
    arrange(sample_id)
  
  rownames(test2) <- test2[,1]  
  
  return(test2)
}


test3 <- run_conversion(cell_settings)