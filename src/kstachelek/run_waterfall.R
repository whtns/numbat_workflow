#!/usr/bin/Rscript

# load required libraries --------------------------------------------------


source("../src/waterfall.R")
source("/home/skevin/storage/TOOLS/waterfall/maxwellbay-waterfall-9d035d9cf75a/src/Waterfall.R")
#source("/home/thor/TOOLS/waterfall/maxwellbay-waterfall-9d035d9cf75a/src/HyperGeoTerms.R")

library(biomaRt)
library(tidyverse)


# load input files --------------------------------------------------------



EXPRESSION_PATH = "../output/FACS_20171031_dshayler_H_sapiens_output/transcripts.tpm_census_matrix.csv" #this points to your gene expression file
EXPRESSION_PATH = c(EXPRESSION_PATH, "../output/FACS_20170407_dshayler_H_sapiens_output/Dominik_stringtie.tpm_census_matrix.csv") #this points to your gene expression file

GROUP_PATH = "../data/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20171031_dshayler_H_sapiens/Sample_number_attributes_103117.csv"
GROUP_PATH = c(GROUP_PATH, "../data/sc_cone_devel/sc_cone_devel_H_sapiens/FACS_20170407_dshayler_H_sapiens/Dominic_Cell_Age_and_Sort_Method_042017.csv")

subset_params <- read.table("../output/dshayler/20170407_20171031_combined_subset_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

jt_complete = map(EXPRESSION_PATH, read.table, header = TRUE, sep="\t") #the full Dataset is jt_complete

jt_complete = map(EXPRESSION_PATH, read.table, header = TRUE, sep="\t") %>% 
  reduce(merge, by=0, all=TRUE) %>% 
  remove_rownames %>% 
  column_to_rownames(var="Row.names")

# names(jt_complete) <- as.numeric(gsub("^.*_([0-9]+).*$", "\\1", EXPRESSION_PATH))
# jt_complete <- readRDS("/home/thor/single_cell_pipeline/results/sunhye/sunhye_census_matrix.rds")

cell_groups = map_dfr(GROUP_PATH, read.table, header = TRUE, sep = ",")
rownames(cell_groups) <- cell_groups[,1]
# names(cell_groups) <- as.numeric(gsub("^.*_([0-9]+).*$", "\\1", EXPRESSION_PATH))

cone_markers <- c("RXRG", "CRX", "ARR3",  "RORB", "GTF2IRD1", "GNAT2", "THRB", "MDM2")


# load requred functions --------------------------------------------------


ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")

#function that confirms genes meet expression conditions
meetsCondition = function(mat, gtet = 0, nCells=0){
  condition_mat = ifelse(mat>=gtet,1,0)
  meets_condition = apply(condition_mat,1,sum) >= nCells
  return(meets_condition)
}

#assign refseq gene symbol from ensembl transcript id
gene_sym_from_trsid <- function(x, query_col, ...){
  
  query_col <- deparse(substitute(query_col))
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = rownames(x), mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  # x <- rownames_to_column(x, "ensembl_transcript_id")
  # output <-  left_join(as.data.frame(biomaRt_result), x, by = "ensembl_transcript_id")
}

#subset expression data for input into pca and hierarchical clustering takes as input
# 1) dataset, 2) treatment groups
# 3) treatment_overlay 4) day comparisons 5) day overlay
subset_jt <- function(color_param = NULL, jt_df, cell_groups_reduced, subset_params = NULL, gray_controls = FALSE){
  #log transform
  jt_df = log2(jt_df+1)
  
  #meets the expression conditions
  #meets_condition = meetsCondition(jt_df,gtet=0.1,nCells = 10)
  #jt_df = as.data.frame(jt_df[meets_condition,])
  
  if (hasArg(subset_params)){
    for (i in seq_along(subset_params)){
      i <- colnames(subset_params)[i]
      cell_groups_reduced <- dplyr::filter(cell_groups_reduced, cell_groups_reduced[[i]] %in% subset_params[[i]])
    }  
  }
  
  jt_df <- jt_df[,cell_groups_reduced$Sample_ID]
  
  clr_names = unique(cell_groups_reduced[[color_param]])
  clrs = RColorBrewer::brewer.pal(length(clr_names), "Set3")
  names(clrs) <- clr_names
  
  
  # if(gray_controls){
  #   cells_used <- cell_groups_reduced$Treatment[rownames(cell_groups_reduced) %in% colnames(dataset_used)]
  #   idx <- grep("shCtrl", cells_used)
  #   day_col_vect[idx] <- "gray"
  # }
  dataset_used <- list("jt_df" = jt_df, "col_vect" = clrs)

}

#quantifies mito read % and generate vector of colors that matches cell id for pca plotting
color_by_var <- function(var_vect){
  var_vect <- as.vector(var_vect)
  colfunc <- colorRampPalette(c("white", "black"))
  var_col_vect <- colfunc(length(var_vect))[rank(var_vect)]
}


#finds genename in ensemble<->gene association tables
getGene = function(gene, association = NULL, ...) {
  
  gene = paste0("^",gene,"$")
  grep_results <- grep(gene, association$hgnc_symbol)
  transcript_ids <- association$ensembl_transcript_id[grep_results]
  #results <- (as.character(association)[lapply(names(t_to_g_association), function(x) which(gene%in%x))])
  return(transcript_ids)
}

color_by_gene_expr <- function(my_df, dimensions, markers, var_col, out_pdf = "../results/sunhye_expr_val_test.pdf", ...){
  
  #  t_to_g_association <- gene_sym_from_trsid(my_df)
  
  coords = Waterfall.PCA(TPM_Data = my_df, ind_sup_vect = ind_sup_vect,
                         plot_title = "cone marker gene expression", twoDcomp = c(4,5),
                         threeD_plot = FALSE, col = var_col, overlay_col = NULL,
                         pdf_out = out_pdf)#change my_df to dataset_used for specific cell group
  coords = coords$pcaN[,1:2]
  pdf(out_pdf, title="cone marker gene expression")
  for (i in markers){
    
    TargetGene <- my_df[getGene(i, t_to_g_association),] #for all cells, multiple genes just rerun from here down
    #TargetGene <- color_data_used[getGene(i),] #for multiple genes just rerun from here down
    TargetGene = 2^TargetGene
    TargetGene=colSums(TargetGene) #to sum
    d<-density(TargetGene)
    plot(d)
    #for summed row
    trans_vect=rep("black", length(TargetGene))
    trans_vect[which(TargetGene>=50)] = "red"
    trans_vect[which(TargetGene<50 & TargetGene>=10)] = "green"
    trans_vect[which(TargetGene<10 & TargetGene>=5)] = "blue"
    trans_vect[which(TargetGene<5 & TargetGene>=2)] = "yellow"
    
    plot(coords,col = trans_vect,main = i, pch=19)
    text(coords[,1],coords[,2],names(TargetGene), cex=0.3)
  }
  dev.off()
}


quantify_mito <- function(dataset_used){
  
  genom_location <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol', 'chromosome_name'),
                          filters = 'ensembl_transcript_id', values = rownames(dataset_used), mart = ensembl)
  
  mito_dataset <- rownames_to_column(dataset_used, "ensembl_transcript_id") %>%
    gather("sample", "read_count", -ensembl_transcript_id) %>%
    inner_join(genom_location, by = "ensembl_transcript_id") %>%
    dplyr::filter(read_count > 0) %>%
    group_by(sample, chromosome_name) %>%
    mutate(expr_trans = sum(read_count)) %>%
    group_by(sample, chromosome_name, expr_trans) %>%
    summarise()
  
  mito_2 <- mito_dataset %>%
    ungroup() %>%
    dplyr::filter(chromosome_name == "MT") %>%
    dplyr::select(-chromosome_name) %>%
    rename(mito_exprs = expr_trans)
  
  mito_3 <- mito_dataset %>%
    ungroup() %>%
    group_by(sample) %>%
    summarise(expr_trans = sum(expr_trans)) %>%
    left_join(mito_2, by = "sample") %>%
    mutate(mito_pct = (mito_exprs / expr_trans)*100) %>%
    arrange(mito_pct) %>%
    .$mito_pct
  
  mito_col_vect <- color_by_var(mito_3)
}

quantify_marker_expr <- function(dataset_used, markers){
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  genom_location <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol', 'chromosome_name'),
                          filters = 'ensembl_transcript_id', values = rownames(dataset_used), mart = ensembl)
  
  marker_dataset <- rownames_to_column(dataset_used, "ensembl_transcript_id") %>%
    gather("sample", "read_count", -ensembl_transcript_id) %>%
    inner_join(genom_location, by = "ensembl_transcript_id")
  
  marker_2 <- marker_dataset %>%
    dplyr::filter(hgnc_symbol %in% markers) %>%
    group_by(sample) %>%
    mutate(marker_expr = sum(read_count))
  
  marker_3 <- marker_dataset %>%
    ungroup() %>%
    group_by(sample) %>%
    summarise(expr_trans = sum(read_count)) %>%
    left_join(marker_2, by = "sample") %>%
    mutate(marker_pct = (marker_expr / expr_trans)*100) %>%
    arrange(marker_pct)
  
  marker_col_vect <- color_by_var(marker_3)
}



# remove bad cells --------------------------------------------------------


#bad cells

bad_cells = bad_cells_martin


# subset jt to remove "bad cells" -----------------------------------------

jt_rm_bad_cells = jt_complete[,-which(colnames(jt_complete)%in%bad_cells)]
dim(jt_complete)
dim(jt_rm_bad_cells)


# #subset jt for specific comparisons -------------------------------------

cycles <- colnames(subset_params)

datasets <- map(cycles, subset_jt, jt_complete, cell_groups, subset_params = subset_params, gray_controls = FALSE)

dataset_used <- subset_jt(color_param = NULL, jt_complete, cell_groups, subset_params = subset_params, 
                                    gray_controls = FALSE)

all_pca <- map(datasets, Waterfall.PCA, plot_title = "PC 3 and 1", twoDcomp = c(3,1),
                        threeD_plot = FALSE, col = dataset_used$col_vect, label = TRUE)


# martin troubleshooting --------------------------------------------------
cells_to_use_733 <- read.table("../../martin_python_scripts/733_cells_to_use.csv")
transcripts_to_use_733 <- read.table("../../martin_python_scripts/733_transcripts_to_use.csv")
cells_to_use_737 <- read.table("../../martin_python_scripts/737_cells_to_use.csv")
transcripts_to_use_737 <- read.table("../../martin_python_scripts/737_transcripts_to_use.csv")

jt_733_test <- jt_complete[which(rownames(jt_complete) %in% transcripts_to_use_733$V2), which(colnames(jt_complete) %in% cells_to_use_733$V2)]
jt_737_test <- jt_complete[which(rownames(jt_complete) %in% transcripts_to_use_737$V2), which(colnames(jt_complete) %in% cells_to_use_737$V2)]
# #get PCA data object ----------------------------------------------------

# for only 2d plots
jt.2PCA = Waterfall.PCA(TPM_Data = jt_complete, plot_title = "PC 3 and 1", twoDcomp = c(3,1),
                        threeD_plot = FALSE, col = dataset_used$col_vect, label = TRUE)


# for 2d and 3d plots
jt.3PCA = Waterfall.PCA(TPM_Data = dataset_used, plot_title = "MT PC3&4&5", twoDcomp = c(3,4),
                        threeDcomp = c(3,4,5), threeD_plot = TRUE, col = day_col_vect, overlay_col = overlay_vect, label = TRUE, callout_cells = NULL)

# for "multiplex" pca plots with 20 dimensions

jt.MPCA = Waterfall.PCA(TPM_Data = dataset_used, ind_sup_vect = ind_sup_vect,
                        plot_title = "multiplex_733_branch_a", twoDcomp = c(1,2),
                        threeD_plot = FALSE, col = day_col_vect, overlay_col = overlay_vect,
                        pdf_out = "../results/sunhye/20170817_pca_1_2_all_cells/only_733_branch_a_color_by_day.pdf",
                        multiplex_plot = TRUE, label = TRUE, callout_cells = NULL)


# subset pca branches -----------------------------------------------------

pc_coords <- as.data.frame(jt.2PCA$pcaN) %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, Dim.1, Dim.2)

branch_pc_coords <- pc_coords %>%
  mutate(branch = ifelse(Dim.1 < 25, "branch_a",
                         ifelse(Dim.2 > 0, "branch_b", "branch_c"))) %>%
  mutate(branch_corrected = ifelse(cell_id %in% c("X62", "X63", "X152"), "branch_a", ifelse(cell_id =="X224", "branch_b", branch)))

branch_cell_groups <- cell_groups %>%
  rownames_to_column("cell_id") %>%
  inner_join(branch_pc_coords) %>%
  dplyr::select(-`Dim.1`,-`Dim.2`)

write.table(branch_cell_groups, "/home/thor/single_cell_pipeline/results/Sunhye_cell_division_day_treatment_and_branch.csv")

keep_cells <- branch_cell_groups %>%
  dplyr::filter(branch_corrected == "branch_a") %>%
  .$cell_id

dataset_used = dataset_used[,which(colnames(dataset_used)%in%keep_cells)]
subset_jt(dataset_used, group_selection = c("sh733", "sh737", "shCtrl"),
          group_overlay   = c("shCtrl"),
          day_selection   = c("Day_4", "Day_6", "Day_8", "Day_12"),
          day_overlay     = c("none"), gray_controls = TRUE)

jt.2PCA_branch = Waterfall.PCA(TPM_Data = dataset_used, ind_sup_vect = ind_sup_vect, plot_title = "MT PC3&4&5", twoDcomp = c(1,2),
                               threeD_plot = FALSE, col = treat_col_vect)

ggplot(branch_pc_coords, aes(x = Dim.1, y = Dim.2, color = branch_corrected)) + geom_point(shape = 32) + geom_text(aes(label = cell_id)) + geom_text(aes(label = ))


# quantify mitochondrial percentage ---------------------------------------

mito_4 <- quantify_mito(dataset_used)

ggplot(mito_4, aes(mito_pct)) + geom_histogram(binwidth = 0.5)

exaxis <- mito_4$sample
ggplot(mito_4, aes(sample)) + geom_bar(aes(weight = mito_pct)) +
  scale_x_discrete(limits = exaxis) +
  theme(axis.text.x = element_text(size=4, angle = 45, vjust = 1, hjust=1)) +
  labs(y= "mitochondrial percentage")

keep_cells <- mito_4$sample

#=======================================

#Hierarchical Clustering
#Use NULL argument for nClusters for unbiased tree cutting otherwise give positive integer. nClusters = NULL by default.
grx = Waterfall.Cluster(dataset_used) #all
grx = Waterfall.Cluster(dataset_used, nClusters = 3,cell.cols = day_col_vect) #all, color by day
grx = Waterfall.Cluster(dataset_used[,cell_groups_reduced[,"Treatment"] != "shCtrl"], nClusters = 3,cell.cols = day_col_vect[cell_groups_reduced[,"Treatment"] != "shCtrl"]) #no control, color by day
grx_wo_branches = Waterfall.Cluster(dataset_used, cell.cols = mito_col_vect)

#Creating color variable from cluster output
jt.cols = grx[,1]

#Removes the side branches by their cluster # (as assigned in Waterfall.Cluster())
#colorRmIdx = which(grx[,2] %in% c(3))

#jt = jt[,-colorRmIdx]
#jt_tmp = jt_tmp[,-colorRmIdx]
#jt.cols = jt.cols[-colorRmIdx]

#Plots in pseudotemporal space
#color plot by day

pdf_733 = "../../martin_python_scripts/733_pca_4_2_transcripts_cells_specified.pdf"
pdf_737 = "../../martin_python_scripts/737_pca_3_1_transcripts_cells_specified.pdf"

jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, ind_sup_vect = ind_sup_vect, col = day_col_vect, overlay_col = overlay_vect, twoDcomp = c(3,1),
                             plot_title = "737_3_1_20171003", nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE,
                             label=TRUE, twoD_plot = TRUE, threeD_plot = FALSE, pdf_out = pdf_733)

#color plot by treatment
jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = treat_col_vect, plot_title = " Run1",
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE,
                             label=TRUE,threeD_plot = FALSE)



#Sorts TPM data by pseudotime and adds dummy variable "Pseudotime" to gene list
jt.cols = jt.cols[rownames(jt.PT)]
jt_tmp = jt_tmp[,rownames(jt.PT)]
jt_tmp["Pseudotime",] = jt.PT[,1]

jt = jt[,rownames(jt.PT)]
jt["Pseudotime",] = jt.PT[,1]

#example of gene-specific expression visualization with HMM raster
#beta is the exponentiation to return to raw FPKM vals (not log)
#gene_txt parameter takes string and prints as title for our graph
#unit the bin size for hmm cmoputation. getGene() takes the gene
PTgeneplot(getGene("SYK",xloctogene),jt,col=jt.cols,gene_txt = "SYK",hmm=TRUE,beta = 10,unit = 3,apply_to_graph=TRUE)

#ranks genes by correlation coefficient with pseudotime
pt_cor = ptCor(jt_tmp,jt_tmp["Pseudotime",],beta = 10,method = "spearman")

#computes an hmm raster across a number of input genes (in this case the top 150 [1:150])
#organizes genes by onset time, showing activation at later and later PT times
rast = cascadingExpression(jt_tmp[c(names(pt_cor$up)[1:150],"Pseudotime"),],
                           beta=10,unit=3,on_prop = .3,avg_mindist_prop = .25,
                           write_img = TRUE,cell_dim = 3.5, asn_table = xloctogene,cex = .65)

