#!/usr/bin/Rscript
#source("/dataVolume/storage/single_cell_pipeline/src/waterfall.R")
source("~/WaterfallDS.R")
source("/dataVolume/storage/TOOLS/waterfall/maxwellbay-waterfall-9d035d9cf75a/utils/baseB.R")
library(biomaRt)
library(dplyr)
suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86


# ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")

#quantile transform
quantTransform = function(mat){
  quant_mat=matrix(0,nrow=nrow(mat), ncol=ncol(mat))
  rownames(quant_mat)=rownames(mat)
  colnames(quant_mat)= colnames(mat)
  for (i in 1:nrow(mat)){
    v = as.numeric(mat[i,])
    q = quantile(v)
    print (q)
    v[v>=q[1] & v<=q[2]] = 1
    v[v>q[2] & v<=q[3]] = 2
    v[v>q[3] & v<=q[4]] = 3
    v[v>q[4] & v<=q[5]] = 4
    quant_mat[i,]= v
  }
  return(quant_mat)
}

#function that relates genomic features for conversion i.e. ensembl to UCSC
# gene_sym_from_trsid <- function(matrix, mart){
#   
#   biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
#                             'ensembl_transcript_id', values = rownames(matrix), mart = mart)
#   biomaRt_result <- as.data.frame(biomaRt_result, stringsAsFactors = FALSE)
#   gene = biomaRt_result$hgnc_symbol
#   ens_trans = biomaRt_result$ensembl_transcript_id
#   names(ens_trans) = gene
#   assign("t_to_g_association",ens_trans,envir = .GlobalEnv)
#   return(ens_trans)
# }

gene_sym_from_trsid <- function(tx_df){
  
  txids <- rownames(tx_df)

  txs <- transcripts(edb, filter = TxIdFilter(txids), columns = c("symbol"))
  
  names(txids) <- txs$symbol
  
  return(txids)
}


print_wilcoxon_rank_sum_table <- function(x, grp_comp){
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = rownames(x), mart = ensembl)
  biomaRt_result <- data.frame(biomaRt_result, stringsAsFactors = FALSE)
  output <- setDT(as.data.frame(x), keep.rownames = TRUE)[]
  output <- rename(output, ensembl_transcript_id = rn)
  output <-  merge(as.data.frame(biomaRt_result), output, by = "ensembl_transcript_id") 
  filename = paste0(grp_comp,"_wilcox_raw_count",".csv")
  write.table(output,filename, sep="\t", quote=FALSE, row.names = FALSE)
}



#function that confirms genes meet expression conditions 
meetsCondition = function(mat, gtet = 0, nCells=0){
  condition_mat = ifelse(mat>=gtet,1,0)
  meets_condition = apply(condition_mat,1,sum) >= nCells
  return(meets_condition)
  
}


#finds genename in ensemble<->gene association tables
getGene = function(gene, association = t_to_g_association) {
  gene = paste0("^",gene,"$")
  grep_results <- grep(gene, names(association))
  transcript_ids <- unname(association[grep_results])
  #results <- (as.character(association)[lapply(names(t_to_g_association), function(x) which(gene%in%x))])
  return(transcript_ids)
}


EXPRESSION_PATH = "/dataVolume/storage/single_cell_pipeline/output/FACS_20171031_dshayler_H_sapiens_output/stringtie_transcripts.tpm_census_matrix.csv" #this points to your gene expression file 
GROUP_PATH = "/dataVolume/storage/single_cell_pipeline/data/sc_cone_devel/sc_cone_devel_H_sapiens/2_seq_dshayler/2_seq_meta_111517.csv"
out_PDF="~/single_cell_pipeline/results/dominic/TwoSeq_no_RodClus_050418.pdf"
pdf(out_PDF, title="TwoSeq_no_RodClus_050418")
cell_groups = read.csv(GROUP_PATH, header = TRUE, stringsAsFactors = FALSE)  
cell_groups = data.frame(lapply(cell_groups, as.character))
rownames(cell_groups) = cell_groups[,1]
#for merged data set
#EXPRESSION_PATH = "/dataVolume/storage/single_cell_pipeline/output/merged_analyses/dshayler_2_seq_w_FACS_20171031_sunlee_day0_census_matrix.csv" #this points to your gene expression file 
#GROUP_PATH = "/dataVolume/storage/single_cell_pipeline/output/merged_analyses/dshayler_2_seq_w_FACS_20171031_sunlee_day0_census_matrix_meta.csv"
#out_PDF="~/single_cell_pipeline/results/dominic/TwoSeq_021218.pdf"
#pdf(out_PDF, title="TwoSeq_021218")
#cell_groups = read.csv(GROUP_PATH, header = TRUE)  
#cell_groups = data.frame(lapply(cell_groups, as.character))
#rownames(cell_groups) = cell_groups[,1]

#poor_reads <- filter(cell_groups, Poor_Read_Number == "Y") %>% 
#  .$row.names
poor_reads <- rownames(cell_groups[cell_groups["Poor_Read_Number"]=="Y" ,])

moderate_alignment <- rownames(cell_groups[cell_groups["Moderate_Alignment"]=="Y" ,])
rod_cells<-rownames(cell_groups[cell_groups["Rod_Cells"]=="Y" ,])
possible_rods<-rownames(cell_groups[cell_groups["Possible_Rods"]=="Y" ,])
c1<- rownames(cell_groups[cell_groups["Collection_Method"]=="C1",])
outliers<- rownames(cell_groups[cell_groups["Outliers"]=="Y",])
non_pr<-rownames(cell_groups[cell_groups["Non_Photoreceptors"]=="Y",])
vsx<-rownames(cell_groups[cell_groups["VSX2_Outlier"]=="Y",])
cluster_rod<-rownames(cell_groups[cell_groups["X9_Cluster_Green_Rods"]=="Y",])

#Loads transcriptome data
jt_complete = read.table(EXPRESSION_PATH,header = TRUE, sep="\t") #the full Dataset is jt_complete
dim(jt_complete)
jt = jt_complete[,-which(colnames(jt_complete) %in% poor_reads)] #change end value for removed cells
jt = jt[,-which(colnames(jt) %in% c1)]
jt = jt[,-which(colnames(jt) %in% outliers)]
jt = jt[,-which(colnames(jt) %in% vsx)]
#jt = jt[,-which(colnames(jt) %in% cluster_rod)]
#jt = jt[,-which(colnames(jt) %in% non_pr)]
dim(jt)
#rownames(jt) = jt[,1]
#jt = jt[,-1]



#log transform
#jt.l=log2(jt_complete+1) #for raw data no removal
jt.l = log2(jt+1) #remove bad cells 

#Get ensemble ID associations for all data (output t_to_g_association)
t_to_g_association <- gene_sym_from_trsid(jt.l)

#meets the expression conditions
jt.l[is.na(jt.l)] <- 0

meets_condition = meetsCondition(jt.l,gtet=0.1,nCells = 10)
jt.l_tmp = as.data.frame(jt.l[meets_condition,])
dim(jt.l_tmp)

cell_groups_reduced = cell_groups[colnames(jt.l_tmp),]

#change cell_groups_reduced
jt.l_without_C1 = jt.l_tmp[,cell_groups_reduced[,"Collection_Method"] != "C1"]
jt.l_without_FACS = jt.l_tmp[,cell_groups_reduced[,"Collection_Method"] != "FACS"]
jt.l_FW15 = jt.l_tmp[,cell_groups_reduced[,"Fetal_Age"] =="15"]
jt.l_retina6=jt.l_tmp[,cell_groups_reduced[,"Collection_Group"] =="6"]
#Change cell_groups_reduced without gene exclusion (For coloring plots)
color_without_C1 = jt.l[,cell_groups_reduced[,"Collection_Method"] != "C1"]
color_without_FACS = jt.l[,cell_groups_reduced[,"Collection_Method"] != "FACS"]
color_FW15 = jt.l[,cell_groups_reduced[,"Fetal_Age"] =="15"]
color_retina6=jt.l[,cell_groups_reduced[,"Collection_Group"] =="6"]

dataset_used = jt.l_tmp
color_data_used= color_without_C1
cell_groups_reduced = cell_groups[colnames(dataset_used),] #groups for raw plot no cells removed

#pull expression for a given gene across all cells (using getGene)
TargetGene <- jt.l_tmp[getGene("NRL"),]
TargetGene = 2^TargetGene

#TargetGene <- quantTransform((TargetGene))

age_cols = c("blue","green","yellow","red","black")
names(age_cols) = c("13","15","16","18","19")
group_cols = c("purple", "green", "blue","yellow","red","pink","orange","cyan","black")
names(group_cols)=c("1","2","3","4","5","6","7","8","9")
seq_cols = c("red","blue")
names(seq_cols)= c("1","2")
#for combined
name_cols = c("red", "blue")
names(name_cols)=c("DS", "SHL")

#for trajectories
trajectory1_cols= c("lightblue","black","lightgreen")
trajectory2_cols=c("lightblue","black","darkgreen","pink")

age_col_vect = age_cols[cell_groups_reduced$Fetal_Age]
group_col_vect = group_cols[cell_groups_reduced$Collection_Group]
seq_col_vect = seq_cols[cell_groups_reduced$Seq_Number]
name_col_vect = name_cols[cell_groups_reduced$collector]

traj1_col_vect <- cell_groups_reduced$Cluster_Color[cell_groups_reduced$Cluster_Color %in% trajectory1_cols]
traj1_col_vect <- gsub("lightgreen", "green", traj1_col_vect)


traj2_col_vect <- cell_groups_reduced$Cluster_Color[cell_groups_reduced$Cluster_Color %in% trajectory2_cols]
traj2_col_vect <- gsub("darkgreen", "green", traj2_col_vect)

#vector list of cells for trajectory
# T1_list<- cell_groups_reduced %>%
#   filter(Cluster_Color %in% c("lightblue", "black", "lightgreen"))

T1_list<- cell_groups_reduced %>%
  dplyr::filter(Cluster_Color %in% c("lightblue", "black", "lightgreen")) %>% 
  pull(Sample_ID) %>% 
  as.character()

# rownames(T1_list)<-T1_list[,1]
# T1_list<-rownames(T1_list)

T2_list<- cell_groups_reduced %>%
  dplyr::filter(Cluster_Color %in% c("lightblue","black","darkgreen","pink")) %>% 
  pull(Sample_ID) %>% 
  as.character()

#color cells of interest in whole plot
col_vect=rep("black",ncol(jt.l_tmp))
names(col_vect)=colnames(jt.l_tmp)

cells_of_interest= which(names(col_vect) %in% rod_cells)
group_2= which(names(col_vect) %in% possible_rods)
col_vect[cells_of_interest]="red"
col_vect[group_2]="blue"
#run Waterfall.PCA() from new waterfall script 05/17/2017

#get PCA data object w/ cool plot!!
jt.PCA = Waterfall.PCA(TPM_Data = jt.l_tmp, col = jt.cols, plot_title = "030618_All_cells", twoDcomp = c(1,2)) 
jt.PCA = Waterfall.PCA(TPM_Data = jt.l_tmp, col = age_col_vect, plot_title = "Ages", twoDcomp = c(3,2)) 
jt.PCA = Waterfall.PCA(TPM_Data = jt.l_tmp, col = group_col_vect, plot_title = "Retina", twoDcomp = c(3,2)) 
jt.PCA = Waterfall.PCA(TPM_Data = jt.l_tmp, col = seq_col_vect, plot_title = "SeqGrp", twoDcomp = c(3,2)) 
#PCA data from isolated cell group 
jt.PCA = Waterfall.PCA(TPM_Data = dataset_used, col = sort_method_vect, plot_title = " Run1", twoDcomp = c(1,2)) 
#get first 5 PCA coordinates for all cells 
pca_coord_5 <- jt.PCA$pcaN
pc_correlation <- jt.PCA$pc.cor
pos_pc_correlation <- pc_correlation$pc2$p
neg_pc_correlation <- pc_correlation$pc2$n 
#positive correlation genes
pos_corr_df <- data.table::setDT(as.data.frame(pos_pc_correlation, stringsAsFactors = FALSE), keep.rownames = TRUE)[]
pos_corr_df <- rename(pos_corr_df, ensembl_transcript_id = rn)
pos_t_to_g_df <- data.frame(gene_id = names(t_to_g_association), ensembl_transcript_id = t_to_g_association)

pos_corr_df <- pos_t_to_g_df %>%
  full_join(pos_corr_df, by = "ensembl_transcript_id") %>%
  arrange(desc(pos_pc_correlation))
write.table(pos_corr_df,"~/single_cell_pipeline/results/dominic/2seq_pc2_RemovedRodNonPR_POS_Waterfallgenelist.csv", sep="\t", quote=FALSE, row.names = FALSE)

#Negative correlation genes
neg_corr_df <- data.table::setDT(as.data.frame(neg_pc_correlation, stringsAsFactors = FALSE), keep.rownames = TRUE)[]
neg_corr_df <- rename(neg_corr_df, ensembl_transcript_id = rn)
neg_t_to_g_df <- data.frame(gene_id = names(t_to_g_association), ensembl_transcript_id = t_to_g_association)

neg_corr_df <- neg_t_to_g_df %>%
  full_join(neg_corr_df, by = "ensembl_transcript_id") %>%
  arrange(neg_pc_correlation)
write.table(neg_corr_df,"~/single_cell_pipeline/results/dominic/2seq_pc2_RemovedRodNonPR_Neg_Waterfallgenelist.csv", sep="\t", quote=FALSE, row.names = FALSE)

#=======================================
#=======================================
#=======================================
#stop here for PCA plot

#Hierarchical Clustering
#Use NULL argument for nClusters for unbiased tree cutting otherwise give positive integer. nClusters = NULL by default. 
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL) #Raw run, nothing removed or colored
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = age_col_vect ) #all, color by fetal age
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = group_col_vect ) #all, color by collection group 
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = seq_col_vect )

#for mixed data set
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = name_col_vect )
#grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = quantile_vect ) # all, color by gene 
#For REMOVING SPECIFIC GROUPS
#grx = Waterfall.Cluster(dataset_used, nClusters = NULL) #Raw run, nothing removed or colored

#grx = Waterfall.Cluster(dataset_used, nClusters = NULL, cell.cols = age_col_vect ) #color by fetal age
#grx = Waterfall.Cluster(dataset_used, nClusters = NULL, cell.cols = group_col_vect ) # color by collection group 


#Creating color variable from cluster output
jt.cols = grx[,1]

#Removes the side branches by their cluster # (as assigned in Waterfall.Cluster())
#colorRmIdx = which(grx[,2] %in% c(3))

#jt = jt[,-colorRmIdx]
#jt_tmp = jt_tmp[,-colorRmIdx]
#jt.cols = jt.cols[-colorRmIdx]

#Plots in pseudotemporal space
jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = jt.cols, plot_title = "TwoSeq_Remove_rodcluster_050418.pdf", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)


jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = age_col_vect, plot_title = "TwoSeq_Remove_rodcluster_050418 by Fetal Age", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = group_col_vect, plot_title = "TwoSeq_Remove_rodcluster_050418 by Retina", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = seq_col_vect, plot_title = "TwoSeq_Remove_rodcluster_050418 by Seq", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)


#for trajectory
jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = traj2_col_vect, plot_title = "TwoSeq_Remove_rodcluster_050418.pdf", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE, threeD_plot = TRUE, Keep=T2_list)

jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = traj1_col_vect, plot_title = "TwoSeq_Remove_rodcluster_050418.pdf", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE, threeD_plot = TRUE, Keep=T1_list)

#for combined data set
jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = name_col_vect, plot_title = "TwoSeq_Remove_VSX_hi_AddSHL_cells_021518 by person", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE)

#for specific cells of interest 
#jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = col_vect, plot_title = "all cells outliers and rods removed colored lower reads 070617", 
                             #nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             #label=TRUE,threeD_plot = FALSE)
#for specific gene, label below
#jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = mid_cells, plot_title = "all cells outliers removed 061917 Sort by RXRG", 
                            # nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                            # label=TRUE,threeD_plot = FALSE)
#WITH CELLS REMOVED
#jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = jt.cols, plot_title = "FACS outliers and rods removed colored lower reads 070617", 
                           #  nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             #label=TRUE,threeD_plot = FALSE)

#jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = sort_method_vect, plot_title = "FACS outliers and rods removed colored lower reads 070617 Sort by C1 vs FACS", 
                            # nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                            # label=TRUE,threeD_plot = FALSE)

#jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = age_col_vect, plot_title = "FACS outliers and rods removed colored lower reads 070617 Sort by Fetal Age", 
                             #nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             #label=TRUE,threeD_plot = FALSE)

#jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = group_col_vect, plot_title = "FACS outliers and rods removed colored lower reads 070617 Sort by Retina", 
                            # nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                            # label=TRUE,threeD_plot = FALSE)

#coloring and plotting gene for each transcript at once 
markers <- c("SKP2", "PROM1", "OTX2","NES", "RXRG", "CRX", "ARR3","RORB", "GTF2IRD1", "MYCN", "NXNL1", "NCL","BSG1", "CDKN1B", "RB1", "OLFM3","GNAT2", "THRB","LHX2", "RLBP1", "CCND1", "MDM2", "SYK", "OPN1SW", "VSX2", "SOX2","CRB1","CRB2","POU4F2", "PAX6", "PROX1", "OPN1MW","OPN1LW", "CRB3","NRL","NR2E3", "PDE6B", "CNGA1","CNGB1", "GNAT1", "GNB1","RHO","PDE6A")
photo_markers <- c("PDE6C", "PDE6H", "PDE6A", "PDE6B", "PDE6G", "CNGA3", "CNGB3", "CNGA1", "CNGB1", "GRK7", "GRK1", "SAG", "GNB3", "GNGT2", "GNB1", "GNGT1", "GUCY2D", "GUCY2F", "GUCA1C", "GUCA1A", "GUCA1B")
rod_markers<-c("NRL","NR2E3", "PDE6B", "CNGA1","CNGB1", "GNAT1", "GNB1","RHO","PDE6A")

#skin_markers<-c("KRT16", "KRT2", "LCE1B", "LCE5A", "SPRR2A", "BNC1", "EMP1", "RBP2", "ABCA12", "CLDN1","SFN")
coords = Waterfall.PCA(jt.l_tmp,twoDcomp = c(1,2))#change jt.l_tmp to dataset_used for specific cell group
coords = coords$pcaN[,1:2]

#FOR SIJIA
sijia_markers<-c("KDM5A","ACO2","BCKDHA","DLAT","DLST","DMPK","GPT2","GRPEL1","GLU1","IDH3B","MAP2K3","MAP2K7","MAPK4","MDH1","MCAT","MLYCD","MTOR","MUL1","NME1","PGD","PIM1","TIMM22","TIMM50","TOMM6","TOMM7","TOMM22","VDAC1","MRPL9","MRPL11","MRPL12","MRPL17","MRPS23","MRPL53","MRPL55","GADD45GIP1","NIT1","PMPCB","TXN2","SOD2","COQ5","COQ7","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA9","NDUFA10","NDUFA13","NDUFB8","NDUFB10","NDUFS2","NDUFS7","NDUFV1","SDHB","SDHC","SDHD","UQCRQ","UQCRC1","UQCRFS1","COA3","COX6A1", "COX6B1", "COX4I1","COX4I2","COX7A2","COX14","COX17","ATP5A1","ATP1A1","ATP1B1","ATP5I","ATP6V1D","ATP6V1F","MFN2","OPA1","PARK7","CHCHD3","DNAJC11","GHITM","MSTO1","REEP1","MRPL58")
out_PDF="/dataVolume/storage/single_cell_pipeline/results/dominic/032118no_vsx_photo_markers.pdf"
pdf(out_PDF, title="032118no_vsx_photo_markers")
color_exp_plot<- plot_cell_colors(sijia_markers, jt)
dev.off()

plot_cell_colors <- function(marker_set, dataset){
  cell_colors <- data.frame("cells" = colnames(dataset))
  for (i in marker_set){
    TargetGene <- dataset[getGene(i),] #for all cells, multiple genes just rerun from here down
    #TargetGene <- color_data_used[getGene(i),] #for multiple genes just rerun from here down
    #TargetGene = (2^TargetGene)-1
    TargetGene=colSums(TargetGene) #to sum
    d<-density(TargetGene)
    plot(d)
    #for summed row
    trans_vect=rep("black", length(TargetGene))
    trans_vect[which(TargetGene>=200)] = "cyan"
    trans_vect[which(TargetGene<200 & TargetGene>=100)] = "red"
    trans_vect[which(TargetGene<100 & TargetGene>=50)] = "orange"
    trans_vect[which(TargetGene<50 & TargetGene>=25)] = "yellow"
    trans_vect[which(TargetGene<25 & TargetGene>=10)] = "green"
    trans_vect[which(TargetGene<10 & TargetGene>=5)] = "blue"
    trans_vect[which(TargetGene<5 & TargetGene>=2)] = "purple"
    trans_vect[which(TargetGene<2 & TargetGene>=0.1)] = "pink"
    plot(coords,col = trans_vect,main = i, pch=19)
    text(coords[,1],coords[,2],names(TargetGene), cex=0.3) 
    cell_colors <- cbind(cell_colors, i = trans_vect)
  }
  
  colnames(cell_colors) <- c("cells", marker_set)
  return(cell_colors)

}

color_exp_plot<- plot_cell_colors(photo_markers, jt)
write.csv(color_exp_plot, file = "/dataVolume/storage/single_cell_pipeline/results/dominic/113017_marker_gene_colors.csv")
dev.off()

  #for individual transcripts
plot_transcript_colors <- function(marker_set, dataset){
  transcript_colors <- data.frame("cells" = colnames(dataset))
  for (i in marker_set){
    TargetGene <- dataset[getGene(i),] #for all cells, multiple genes just rerun from here down
    #TargetGene <- color_data_used[getGene(i),] #for multiple genes just rerun from here down
    #TargetGene = (2^TargetGene)-1
    d<-density(TargetGene)
    plot(d)
    #for summed row
    trans_vect=rep("black", length(TargetGene))
    trans_vect[which(TargetGene>=200)] = "cyan"
    trans_vect[which(TargetGene<200 & TargetGene>=100)] = "red"
    trans_vect[which(TargetGene<100 & TargetGene>=50)] = "orange"
    trans_vect[which(TargetGene<50 & TargetGene>=25)] = "yellow"
    trans_vect[which(TargetGene<25 & TargetGene>=10)] = "green"
    trans_vect[which(TargetGene<10 & TargetGene>=5)] = "blue"
    trans_vect[which(TargetGene<5 & TargetGene>=2)] = "purple"
    trans_vect[which(TargetGene<2 & TargetGene>=0.1)] = "pink"
    plot(coords,col = trans_vect,main = i, pch=19)
    text(coords[,1],coords[,2],names(TargetGene), cex=0.3) 
    cell_colors <- cbind(cell_colors, i = trans_vect)
  }
  colnames(transcript_colors) <- c("cells", marker_set)
  return(transcript_colors)
}

trans_exp_plot<- plot_transcript_colors(markers, jt)
for (i in 1:nrow(TargetGene)){  
  exp = TargetGene[i,]
  trans_vect = rep("black",length(exp))
  trans_vect[which(TargetGene>=200)] = "cyan"
  trans_vect[which(TargetGene<200 & TargetGene>=100)] = "red"
  trans_vect[which(TargetGene<100 & TargetGene>=50)] = "orange"
  trans_vect[which(TargetGene<50 & TargetGene>=25)] = "yellow"
  trans_vect[which(TargetGene<25 & TargetGene>=10)] = "green"
  trans_vect[which(TargetGene<10 & TargetGene>=5)] = "blue"
  trans_vect[which(TargetGene<5 & TargetGene>=2)] = "purple"
  trans_vect[which(TargetGene<2 & TargetGene>=0.1)] = "pink"
  nm = rownames(TargetGene)[i]
  plot(coords,col = trans_vect,main = nm, pch=19)
  text(coords[,1],coords[,2],colnames(TargetGene), cex=0.3) 
}



#Sorts TPM data by pseudotime and adds dummy variable "Pseudotime" to gene list
jt.cols = jt.cols[rownames(jt.PT)]
jt.lt = jt.l[,rownames(jt.PT)]
jt.lt["Pseudotime",] = jt.PT[,1]

jt = dataset_used[,rownames(jt.PT)]
jt["Pseudotime",] = jt.PT[,1]


out_PDF="/dataVolume/storage/single_cell_pipeline/results/dominic/110917_remove_poor_reads_GENES.pdf"
pdf( out_PDF, title="110917_remove_poor_reads_GENES")
#example of gene-specific expression visualization with HMM raster
#beta is the exponentiation to return to raw FPKM vals (not log)
#gene_txt parameter takes string and prints as title for our graph
#unit the bin size for hmm cmoputation. getGene() takes the gene 
PTgeneplot(getGene("SOX2"),jt,col=jt.cols,gene_txt = "SOX2",hmm=TRUE,beta = 2,unit = 3,apply_to_graph=TRUE)


#ranks genes by correlation coefficient with pseudotime
#pt_cor = ptCor(jt.l,jt.l_tmp["Pseudotime",],beta = 10,method = "spearman")
pt_cor = ptCor(jt.lt,jt.lt["Pseudotime",],beta = 2,method = "spearman")
pt_cor_up <- data.frame(ensembl_id =names(pt_cor$up), correlation=pt_cor$up)
pt_cor_dn <- data.frame(ensembl_id =names(pt_cor$dn), correlation=pt_cor$dn)
write.table(as.data.frame(pt_cor_up),"./FACS_cell_psuedo_correl.csv", sep = ",")
write.table(as.data.frame(pt_cor_dn), "./FACS_cell_psuedo_correl.csv", sep = ",", append = TRUE)
test <-  read.table("./FACS_cell_psuedo_correl.csv", sep = ",")
#computes an hmm raster across a number of input genes (in this case the top 150 [1:150])
#organizes genes by onset time, showing activation at later and later PT times 
#rast = cascadingExpression(jt.l_tmp[c(names(pt_cor$up)[1:150],"Pseudotime"),],
                           #beta=10,unit=3,on_prop = .3,avg_mindist_prop = .25,
                           #write_img = TRUE,cell_dim = 3.5, asn_table = t_to_g_association,cex = .65)
rast = cascadingExpression(jt.lt[c(names(pt_cor$up)[1:150],"Pseudotime"),],
                           beta=10,unit=3,on_prop = .3,avg_mindist_prop = .25,
                           write_img = TRUE,cell_dim = 3.5, asn_table = t_to_g_association,cex = .65)
rast = cascadingExpression(jt.lt[c(names(pt_cor$dn)[1:150],"Pseudotime"),],
                           beta=10,unit=3,on_prop = .3,avg_mindist_prop = .25,
                           write_img = TRUE,cell_dim = 3.5, asn_table = t_to_g_association,cex = .65)


# look up number of cells expressing a given gene -------------------------

audit_gene_counts <- function(gene, g_c_matrix, thresh){
  
  trs_vec <- getGene(gene)
  
  lapply(trs_vec, function(x){g_c_matlsrix[x,]>thresh})
  
}


match 

  