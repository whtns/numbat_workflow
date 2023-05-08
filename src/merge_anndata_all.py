#!/usr/bin/env python 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import patchworklib as pw
import scanpy as sc
import scvi
from pathlib import Path
import anndata
import pandas as pd
from matplotlib.pyplot import rc_context
import os 
import glob
import numpy as np
from scvi.model.utils import mde
import pymde
import re
from plotnine import ggplot, geom_bar, aes, coord_flip

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)
  
def plot_my_embedding(adata_input, color, basis="X_mde"):
  myadata = adata_input.copy()
  sc.pl.embedding(
    myadata,
    basis=basis,
    color=color,
    frameon=False,
    ncols=1,
    save = f'_{color}.pdf'
)


adata_all = sc.read_h5ad("output/scanpy/scvi_all.h5ad")

# mitochondrial genes
adata_all.var['mt'] = adata_all.var_names.str.startswith('MT-') 
# ribosomal genes
adata_all.var['ribo'] = adata_all.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
adata_all.var['hb'] = adata_all.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(adata_all, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata_all, min_genes=200)
sc.pp.filter_genes(adata_all, min_cells=3)
# filter for percent mito
adata_all = adata_all[adata_all.obs['pct_counts_mt'] < 20, :]
# filter for percent ribo > 0.05
adata_all = adata_all[adata_all.obs['pct_counts_ribo'] > 5, :]
malat1 = adata_all.var_names.str.startswith('MALAT1')
# we need to redefine the mito_genes since they were first 
# calculated on the full object before removing low expressed genes.
mito_genes = adata_all.var_names.str.startswith('MT-')
hb_genes = adata_all.var_names.str.contains('^HB[^(P)]')

remove = np.add(mito_genes, malat1)
remove = np.add(remove, hb_genes)
keep = np.invert(remove)

adata_all = adata_all[:,keep]

adata_all.obs.cluster_short.replace("ARLIP1", "ARL6IP1")
adata_all = adata_all[~adata_all.obs.clone_opt.isna()]


sc.pp.log1p(adata_all)
sc.tl.leiden(adata_all, resolution = 0.2)
adata_all.obsm["X_mde"] = mde(adata_all.obsm["X_scVI"])

leiden_membership = adata_all.obs.leiden.value_counts()

too_small_leiden_clusters = pd.Series(np.where(leiden_membership < 10))

adata_all = adata_all[~adata_all.obs.leiden.isin(["13", "14", "15", "16", "17"]), :]

plot_my_embedding(adata_all, "leiden")

sc.tl.rank_genes_groups(adata_all, 'leiden', method='wilcoxon')

sc.pl.rank_genes_groups(adata_all, n_genes=10, sharey=False, save = "_leiden_marker_genes.pdf")

sc.tl.rank_genes_groups(adata_all, 'cluster_short', method='t-test')
sc.pl.rank_genes_groups(adata_all, n_genes=1, sharey=False, save = "_short_marker_genes.pdf")

adata_all.obs['cluster_short'] = adata_all.obs.cluster.str.replace(".*_", "")

adata_all.obs['clone_opt'] = adata_all.obs.clone_opt.astype('category')

# write obs data to file for plotting R ggplot2 ------------------------------
adata_all.obs.to_csv("results/adata_all_metadata.csv")

adata_all.write_h5ad("output/scanpy/adata_all_processd.h5ad")

adata_all =  sc.read_h5ad("output/scanpy/adata_all_processd.h5ad")
adata_all.uns['log1p']["base"] = None
# 
# 
# top_genes = list()
# for i in range(0,13):
#   top_genes = top_genes + adata_all.uns['rank_genes_groups']['names'].field(i)[0:1].tolist()
# 
# interesting_genes  = ["RXRG", "DEK", "KIF14", "SOX4", "NEK2"]
# 
# top_genes = top_genes  + interesting_genes
# 
# sc.pl.matrixplot(
#     adata_all,
#     var_names = top_genes,
#     groupby ='leiden',
#     dendrogram = True,
#     log = True,
#     save = "matrixplot.pdf"
# )

plot_my_embedding(adata_all, "sample_id")
plot_my_embedding(adata_all, "leiden")
plot_my_embedding(adata_all, "cluster_short")
plot_my_embedding(adata_all, "clone_opt")

# filter bad cells ------------------------------

cluster_labels = {
  '0' :  "TFF1",
  '1' :  "injured rb",
  '2' :  "G2M",
  '3' :  "injured rod",
  '4' :  "rod",
  '5' :  "immune",
  '6' :  "immune",
  '7' :  "calm1",
  '8' :  "injured cone",
  '9' :  "rpe",
  '10' : "s",
  '11' : "cone",
  '12' : "apoe"
  }
  
adata_all.obs['annotated_leiden'] = adata_all.obs['leiden'].map(cluster_labels)

plot_my_embedding(adata_all, "annotated_leiden")

bad_clusters = ["injured rod", "rod", "immune", "rpe", "apoe"]

adata_bad = adata_all[(adata_all.obs.annotated_leiden.isin(bad_clusters)), :].copy()

plot_my_embedding(adata_bad, "clone_opt")
plot_my_embedding(adata_bad, "sample_id")

adata_filtered = adata_all[~(adata_all.obs.annotated_leiden.isin(bad_clusters)), :].copy()

plot_my_embedding(adata_filtered, "annotated_leiden")
plot_my_embedding(adata_filtered, "cluster_short")

sc.tl.leiden(adata_filtered, resolution = 0.3)

plot_my_embedding(adata_filtered, "leiden")
sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_filtered, n_genes=25, sharey=False, save = "_leiden_marker_genes.pdf")

filtered_cluster_labels = {
  '0' :  "TFF1",
  '1' :  "G2M",
  '2' :  "HSP",
  '3' :  "XIST",
  '4' :  "injured rod",
  '5' :  "bipolar",
  '6' :  "ARL6IP1",
  '7' :  "MT 1",
  '8' :  "unknown 1",
  '9' :  "MT 2",
  '10' : "unknown 2",
  '11' : "cone",
  '12' : "unknown 3"
  }


adata_filtered.obs['annotated_leiden'] = adata_filtered.obs['leiden'].map(filtered_cluster_labels)

plot_my_embedding(adata_filtered, "annotated_leiden")

adata_filtered.write_h5ad("output/scanpy/adata_all_filtered.h5ad")

adata_filtered = sc.read_h5ad("output/scanpy/adata_all_filtered.h5ad")
adata_filtered.uns['log1p']["base"] = None

adata_filtered.obs['cluster_short'] = adata_filtered.obs['cluster_short'].replace("ARL1IP1", "ARL6IP1")

# asdf

leiden_membership = adata_filtered.obs.leiden.value_counts()

cells_w_enough_leiden = adata_filtered.obs.groupby("leiden").filter(lambda x: len(x) > 10)

adata_filtered[adata_filtered.obs.index.isin(cells_w_enough_leiden.index), :]

# remove bad clusters ------------------------------
cluster_dictionary = pd.read_csv("data/cluster_dictionary.csv")

cluster_dictionary['cluster'] = cluster_dictionary["sample_id"] + "_" + cluster_dictionary["abbreviation"]

cluster_dictionary = cluster_dictionary.loc[cluster_dictionary['remove'] == 1]

adata_filtered = adata_filtered[~adata_filtered.obs.cluster.isin(cluster_dictionary.cluster), :]

adata_filtered.obs.to_csv("results/adata_filtered_metadata.csv")

def plot_dend(adata):
  adata_filtered.obs['cluster_short'] = adata_filtered.obs['cluster_short'].replace("ARL1IP1", "ARL6IP1")
  adata_filtered.obs['cluster_short'] = adata_filtered.obs.cluster_short.astype('category')
  sc.tl.dendrogram(adata_filtered, groupby = "cluster_short")
  dend = sc.pl.dendrogram(adata_filtered, groupby = "cluster_short", orientation = "left")
  dend.figure.set_size_inches(18, 13)
  dend.figure.savefig(f'_dend.pdf', dpi=200)

# cell cycle phase ------------------------------
cl = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/"
cl += "180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
cc_prots = pd.read_csv(cl, header=None).squeeze()
s_genes = cc_prots[:43].tolist()
g2m_genes = cc_prots[43:].tolist()


sc.tl.score_genes_cell_cycle(adata_filtered, s_genes = s_genes, g2m_genes = g2m_genes)

plot_my_embedding(adata_filtered, "phase")

plot_my_embedding(adata_filtered, "clone_opt")

sc.tl.leiden(adata_filtered, resolution = 0.2)

leiden_membership = adata_filtered.obs.leiden.value_counts()

cells_w_enough_leiden = adata_filtered.obs.groupby("leiden").filter(lambda x: len(x) > 10)

adata_filtered = adata_filtered[adata_filtered.obs.index.isin(cells_w_enough_leiden.index), :]
plot_my_embedding(adata_filtered, "leiden")

sc.tl.rank_genes_groups(adata_filtered, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_filtered, n_genes=25, sharey=False, save = "_leiden_marker_genes.pdf")


plot_my_embedding(adata_filtered, "cluster_short")

plot_my_embedding(adata_filtered, "sample_id")


def scanpy_plot_each_cluster(input_adata, meta):
  myadata = input_adata.copy()
  # myadata = myadata[~myadata.obs.clone_opt.isna()]
  meta_vals = myadata.obs[meta].unique().tolist()
  # meta_vals.sort()
  
  myadata_meta = myadata.obs[meta]
  print(myadata_meta)

  # meta_vals = natural_sort(meta_vals)
  
  meta_vals = [x for x in meta_vals if str(x) != 'nan']
  print(meta_vals)
  
  for i in meta_vals:
    
    spec_meta = np.where(myadata_meta == i, i, pd.NA)
    
    myadata.obs[meta] = spec_meta
    
    sc.pl.embedding(
        myadata,
        basis="X_mde",
        color=meta,
        frameon=False,
        ncols=1,
        save = f'_{meta}_{i}.png',
        use_raw = False,
        na_in_legend = False,
        na_color = 'gray'
    )

scanpy_plot_each_cluster(adata_all, "leiden")

scanpy_plot_each_cluster(adata_all, "cluster_short")



sc.pl.embedding(
    adata_all,
    basis="X_mde",
    color=top_genes,
    frameon=False,
    ncols=1,
    save = "_top_markers.png", 
    use_raw = False
)
# 
# (ggplot(adata_all.obs, aes('wt', 'mpg', color='factor(gear)'))
#  + geom_point()
#  + stat_smooth(method='lm')
#  + facet_wrap('~gear'))
#  
#  ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +
#   geom_col(position = "fill")

# ggplot(adata_all.obs, aes(x = .data[[clusters]], fill = clone_opt)) +
#     geom_bar(position = "fill") +
#     coord_flip()

(ggplot(adata_all.obs, aes('factor(cluster_short)', fill = 'factor(leiden)'))
    + geom_bar(position = "fill"))
    # + coord_flip())
    
adata_all.obs.plot.bar(x='leiden', stacked=True)


# end ------------------------------


# leiden_clusters = adata_all.obs.leiden.unique().tolist()
# leiden_clusters.sort()
# 
# os.makedirs("figures/violin_adata_all")
# os.makedirs("figures/rank_genes_groups_clone_opt_adata_all")
# 
# for i in leiden_clusters:
#   minor_adata = adata_all[adata_all.obs.leiden == i,:]
#   sc.tl.rank_genes_groups(minor_adata, 'clone_opt', method='t-test')
#   sc.pl.rank_genes_groups(minor_adata, n_genes=10, sharey=False, save = f'_adata_all/{i}_marker_genes.pdf')
#   top_genes = list()
#   
#   n_clones = len(minor_adata.obs.clone_opt.unique())
#   
#   print(n_clones)
#   
#   for j in range(0,n_clones):
#     top_genes = top_genes + minor_adata.uns['rank_genes_groups']['names'].field(j)[0:1].tolist()
# 
#   print(top_genes)
#   sc.pl.violin(
#     minor_adata,
#     keys=top_genes,
#     groupby='clone_opt',
#     rotation=90,
#     log = False,
#     use_raw = True,
#     save = f'_adata_all/{i}.pdf')
