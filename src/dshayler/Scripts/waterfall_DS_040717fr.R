#!/usr/bin/Rscript
source("/home/thor/single_cell_pipeline/src/waterfall.R")
source("/home/thor/TOOLS/waterfall/maxwellbay-waterfall-9d035d9cf75a/utils/baseB.R")
library(biomaRt)
library(dplyr)

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
gene_sym_from_trsid <- function(matrix){
  ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
  biomaRt_result <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'hgnc_symbol'), filters =
                            'ensembl_transcript_id', values = rownames(matrix), mart = ensembl)
  biomaRt_result <- as.data.frame(biomaRt_result, stringsAsFactors = FALSE)
  gene = biomaRt_result$hgnc_symbol
  ens_trans = biomaRt_result$ensembl_transcript_id
  names(ens_trans) = gene
  assign("t_to_g_association",ens_trans,envir = .GlobalEnv)
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
getGene = function(gene,association = t_to_g_association) {
  gene = paste0("^",gene,"$")
  grep_results <- grep(gene, names(association))
  transcript_ids <- unname(association[grep_results])
  #results <- (as.character(association)[lapply(names(t_to_g_association), function(x) which(gene%in%x))])
  return(transcript_ids)
}

#saved Objects
#jt.l
jt.l<-readRDS("/media/thor/storage/single_cell_pipeline/Dominic_Waterfall_Analysis/Saved_R_objects/jt.l")
#jt.l_tmp
jt.l_tmp<-readRDS("/media/thor/storage/single_cell_pipeline/Dominic_Waterfall_Analysis/Saved_R_objects/jt.l_tmp")
#t_to_g_association
t_to_g_association<-readRDS("/media/thor/storage/single_cell_pipeline/Dominic_Waterfall_Analysis/Saved_R_objects/t_to_g_association")
#cell_groups_reduced
cell_groups_reduced<-readRDS("/media/thor/storage/single_cell_pipeline/Dominic_Waterfall_Analysis/Saved_R_objects/cell_groups_reduced")

EXPRESSION_PATH = "/home/thor/single_cell_pipeline/results/dominic_census_matrix.csv" #this points to your gene expression file 
#EXPRESSION_PATH = "/home/thor/single_cell_pipeline/results/Dominik_stringtie.tpm.csv"  #TPM file for census
GROUP_PATH = "/home/thor/single_cell_pipeline/results/Dominic_Cell_Age_and_Sort_Method_042017.csv"
#CONVERSION_FILE_PATH = "E:/glps01/GLPS01/David_Cobrinik/113015_DS_HS_mouse_human_rnaseq/human/cuffnorm/genes.attr_table" #this path will need to be changed to point to your genes.attr_table file 
out_PDF="/home/thor/single_cell_pipeline/Dominic_Waterfall_Analysis/FACS outliers and rods removed colored lower reads 070617.pdf"
pdf( out_PDF, title="FACS outliers and rods removed colored lower reads 070617.pdf")
#bad cells
bad_cells = c("X422","X503", "X468", "X399", "X392", "X398", "X396", "X388", "X382", "X381", "X386", "X387", "X395", "X394","X449", "X499", "X383", "X488", "X486", "X391")
facs_outliers=c(bad_cells,"X397", "X390", "X393")   #22
rod_cells=(c(facs_outliers,"X361","X389","X402","X323","X333","X336","X337","X345","X363","X401","X410","X428","X322","X325","X344","X362","X400","X430","X329","X330","X343","X419","X431","X434","X510"))#46
mid_aligned=c(rod_cells,"X490","X453","X497","X474","X493","X475","X484","X476","X478","X479","X482","X480","X509")
#Loads transcriptome data
jt_complete = read.table(EXPRESSION_PATH,header = TRUE, sep=",") #the full Dataset is jt_complete
dim(jt_complete)
jt = jt_complete[,-which(colnames(jt_complete)%in%rod_cells)] #change end value for removed cells
dim(jt)
rownames(jt) = jt[,1]
jt = jt[,-1]



#log transform
#jt.l=log2(jt_complete+1) #for raw data no removal
jt.l = log2(jt+1) #remove bad cells 

#Get ensemble ID associations for all data (output t_to_g_association)
gene_sym_from_trsid(jt.l)

#meets the expression conditions
meets_condition = meetsCondition(jt.l,gtet=0.1,nCells = 10)
jt.l_tmp = as.data.frame(jt.l[meets_condition,])
dim(jt.l_tmp)



cell_groups = read.csv(GROUP_PATH, header = TRUE, row.names=1, colClasses=c("character","character","character"))
cell_groups_reduced = cell_groups[colnames(jt.l_tmp),]
#change cell_groups_reduced
jt.l_without_C1 = jt.l_tmp[,cell_groups_reduced[,"Sort_Method"] != "C1"]
jt.l_without_FACS = jt.l_tmp[,cell_groups_reduced[,"Sort_Method"] != "FACS"]
jt.l_FW15 = jt.l_tmp[,cell_groups_reduced[,"Fetal_Age"] =="15"]
jt.l_retina6=jt.l_tmp[,cell_groups_reduced[,"Collection_Group"] =="6"]
#Change cell_groups_reduced without gene exclusion (For coloring plots)
color_without_C1 = jt.l[,cell_groups_reduced[,"Sort_Method"] != "C1"]
color_without_FACS = jt.l[,cell_groups_reduced[,"Sort_Method"] != "FACS"]
color_FW15 = jt.l[,cell_groups_reduced[,"Fetal_Age"] =="15"]
color_retina6=jt.l[,cell_groups_reduced[,"Collection_Group"] =="6"]

dataset_used = jt.l_without_C1
color_data_used= color_without_C1
cell_groups_reduced = cell_groups[colnames(dataset_used),] #groups for raw plot no cells removed

#pull expression for a given gene across all cells (using getGene)
TargetGene <- jt.l_tmp[getGene("NRL"),]
TargetGene = 2^TargetGene

#TargetGene <- quantTransform((TargetGene))

sort_cols = c("blue","red")
names(sort_cols) = c("C1","FACS")
age_cols = c("green","orange","black")
names(age_cols) = c("15","16","18")
group_cols = c("purple", "green", "blue","yellow","red","black")
names(group_cols)=c("1","2","3","4","5","6")
mid_cells=c("X490","X453","X497","X474","X493","X475","X484","X476","X478","X479","X482","X480","X509")
#quant_cols=c("black","green","yellow","red")
#quantile_vect=quant_Retina 6 Only outliers and rods removed 060517cols[TargetGene[1,]]
sort_method_vect = sort_cols[cell_groups_reduced$Sort_Method]
age_col_vect = age_cols[cell_groups_reduced$Fetal_Age]
group_col_vect = group_cols[cell_groups_reduced$Collection_Group]

#color cells of interest in whole plot
col_vect=rep("black",ncol(jt.l_tmp))
cells_of_interest=which(colnames(jt.l_tmp)%in% c("X490","X453","X497","X474","X493","X475","X484","X476","X478","X479","X482","X480","X509"))
col_vect[cells_of_interest]="red"
#run Waterfall.PCA() from new waterfall script 05/17/2017

#get PCA data object w/ cool plot!!
#jt.PCA = Waterfall.PCA(TPM_Data = jt.l_tmp, col = sort_method_vect, plot_title = " Run1", twoDcomp = c(1,2)) 
#PCA data from isolated cell group 
jt.PCA = Waterfall.PCA(TPM_Data = dataset_used, col = sort_method_vect, plot_title = " Run1", twoDcomp = c(1,2)) 
#get first 5 PCA coordinates for all cells 
pca_coord_5 <- jt.PCA$pcaN
pc_correlation <- jt.PCA$pc.cor
pos_pc_correlation <- pc_correlation$pc2$p
neg_pc_correlation <- pc_correlation$pc2$n 
#positive correlation genes
pos_corr_df <- setDT(as.data.frame(pos_pc_correlation, stringsAsFactors = FALSE), keep.rownames = TRUE)[]
pos_corr_df <- rename(pos_corr_df, ensembl_transcript_id = rn)
pos_t_to_g_df <- data.frame(gene_id = names(t_to_g_association), ensembl_transcript_id = t_to_g_association)

pos_corr_df <- pos_t_to_g_df %>%
  full_join(pos_corr_df, by = "ensembl_transcript_id") %>%
  arrange(desc(pos_pc_correlation))
write.table(pos_corr_df,"Retina6_pc2_POS_Waterfallgenelist.csv", sep="\t", quote=FALSE, row.names = FALSE)

#Negative correlation genes
neg_corr_df <- setDT(as.data.frame(neg_pc_correlation, stringsAsFactors = FALSE), keep.rownames = TRUE)[]
neg_corr_df <- rename(neg_corr_df, ensembl_transcript_id = rn)
neg_t_to_g_df <- data.frame(gene_id = names(t_to_g_association), ensembl_transcript_id = t_to_g_association)

neg_corr_df <- neg_t_to_g_df %>%
  full_join(neg_corr_df, by = "ensembl_transcript_id") %>%
  arrange(neg_pc_correlation)
write.table(neg_corr_df,"Retina6_pc2_NEG_Waterfallgenelist.csv", sep="\t", quote=FALSE, row.names = FALSE)

#=======================================
#=======================================
#=======================================
#stop here for PCA plot

#Hierarchical Clustering
#Use NULL argument for nClusters for unbiased tree cutting otherwise give positive integer. nClusters = NULL by default. 
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL) #Raw run, nothing removed or colored
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = sort_method_vect) #all, color by sort method
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = age_col_vect ) #all, color by fetal age
grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = group_col_vect ) #all, color by collection group 
#grx = Waterfall.Cluster(jt.l_tmp, nClusters = NULL, cell.cols = quantile_vect ) # all, color by gene 
#For REMOVING SPECIFIC GROUPS
grx = Waterfall.Cluster(dataset_used, nClusters = NULL) #Raw run, nothing removed or colored
grx = Waterfall.Cluster(dataset_used, nClusters = NULL, cell.cols = sort_method_vect) # color by sort method
grx = Waterfall.Cluster(dataset_used, nClusters = NULL, cell.cols = age_col_vect ) #color by fetal age
grx = Waterfall.Cluster(dataset_used, nClusters = NULL, cell.cols = group_col_vect ) # color by collection group 


#Creating color variable from cluster output
jt.cols = grx[,1]

#Removes the side branches by their cluster # (as assigned in Waterfall.Cluster())
#colorRmIdx = which(grx[,2] %in% c(3))

#jt = jt[,-colorRmIdx]
#jt_tmp = jt_tmp[,-colorRmIdx]
#jt.cols = jt.cols[-colorRmIdx]

#Plots in pseudotemporal space
jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = jt.cols, plot_title = "all cells outliers removed 061917.pdf", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = sort_method_vect, plot_title = "all cells outliers removed 061917 Sort by C1 vs FACS", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = age_col_vect, plot_title = "all cells outliers removed 061917 Sort by Fetal Age", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = group_col_vect, plot_title = "all cells outliers removed 061917 Sort by Retina", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)
#for specific cells of interest 
jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = col_vect, plot_title = "all cells outliers and rods removed colored lower reads 070617", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)
#for specific gene, label below
jt.PT = Waterfall.Pseudotime(jt.l_tmp, angle=0, col = mid_cells, plot_title = "all cells outliers removed 061917 Sort by RXRG", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)
#WITH CELLS REMOVED
jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = jt.cols, plot_title = "FACS outliers and rods removed colored lower reads 070617", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = sort_method_vect, plot_title = "FACS outliers and rods removed colored lower reads 070617 Sort by C1 vs FACS", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = age_col_vect, plot_title = "FACS outliers and rods removed colored lower reads 070617 Sort by Fetal Age", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

jt.PT = Waterfall.Pseudotime(dataset_used, angle=0, col = group_col_vect, plot_title = "FACS outliers and rods removed colored lower reads 070617 Sort by Retina", 
                             nclust = 5, seed = 5,scalePCA = F,mst=FALSE,invert_pt = FALSE, 
                             label=TRUE,threeD_plot = FALSE)

#coloring and plotting gene for each transcript at once 
markers <- c("RXRG", "CRX", "ARR3", "POU4F2", "PROX1", "RHO", "RORB", "GTF2IRD1", "NR2E3", "MYCN", "CDKN1B", "RB1", "OLFM3", "GNAT1", "NRL", "GNAT2", "THRB", "LHX2", "CCND1", "MDM2", "SYK", "OPN1SW","OPN1MW","OPN1LW", "VSX2", "SOX2","CRB1","CRB2","CRB3")
rod_markers<-c("NRL","NR2E3", "RHO", "PDE6A", "PDE6B", "CNGA1","CNGB1", "GNAT1", "GNB1")

coords = Waterfall.PCA(jt.l_tmp)#change jt.l_tmp to dataset_used for specific cell group
coords = coords$pcaN[,1:2]
out_PDF="/home/thor/single_cell_pipeline/Dominic Waterfall Analysis/all cells outliers removed 061917.pdf"
pdf( out_PDF, title="all cells outliers removed 061917")


plot_cell_colors <- function(marker_set, dataset){
  cell_colors <- data.frame("cells" = colnames(dataset))
  for (i in marker_set){
    TargetGene <- dataset[getGene(i),] #for all cells, multiple genes just rerun from here down
    #TargetGene <- color_data_used[getGene(i),] #for multiple genes just rerun from here down
    TargetGene = (2^TargetGene)-1
    TargetGene=colSums(TargetGene) #to sum
    d<-density(TargetGene)
    plot(d)
    #for summed row
    trans_vect=rep("black", length(TargetGene))
    trans_vect[which(TargetGene>=100)] = "red"
    trans_vect[which(TargetGene<100 & TargetGene>=50)] = "orange"
    trans_vect[which(TargetGene<50 & TargetGene>=25)] = "yellow"
    trans_vect[which(TargetGene<25 & TargetGene>=10)] = "green"
    trans_vect[which(TargetGene<10 & TargetGene>=5)] = "blue"
    trans_vect[which(TargetGene<5 & TargetGene>=2)] = "purple"
    plot(coords,col = trans_vect,main = i, pch=19)
    text(coords[,1],coords[,2],names(TargetGene), cex=0.3) 
    cell_colors <- cbind(cell_colors, i = trans_vect)
  }
  colnames(cell_colors) <- c("cells", marker_set)
  return(cell_colors)
  
}
color_exp_plot<- plot_cell_colors(markers, jt.l)

dev.off()
  #for individual transcripts
for (i in 1:nrow(TargetGene)){  
  exp = TargetGene[i,]
  trans_vect = rep("black",length(exp))
  
  trans_vect[which(exp>=50)] = "red"
  trans_vect[which(exp<50 & exp>=10)] = "green"
  trans_vect[which(exp<10 & exp>=5)] = "blue"
  
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


out_PDF="/home/thor/single_cell_pipeline/Dominic Waterfall Analysis/GENE PLOTS Census all cells 051917.pdf"
pdf( out_PDF, title="Census all cells GENE PLOTS 051917")
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



# compute wilcoxon rank-sum test for counts -------------------------------
dataset_used = jt.l_tmp
signature_vector <- rownames(cell_groups_reduced)
names(signature_vector) <- cell_groups_reduced$Sort_Method
signature_vector <- signature_vector[names(signature_vector) == "FACS"]
###
wilcox_sig_matrix <- signatureAll(dataset_used, signature_vector)
###
wilcox_adjust <- p.adjust(wilcox_sig_matrix$pval, "BH")
wilcox_2 <- wilcox_sig_matrix %>%
  add_column(p.adjust = wilcox_adjust, .after = 2) %>%
  filter(meanA<Inf, log2FC<Inf, !is.na(log2FC)) 
wilcox_2 <- as.data.frame(wilcox_2)
rownames(wilcox_2) <- wilcox_2$name
print_wilcoxon_rank_sum_table(wilcox_2, "FACS_over_C1")
