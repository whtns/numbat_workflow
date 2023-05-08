#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

def plot_my_embedding(adata_input, color, label):
  myadata = adata_input.copy()
  sc.pl.embedding(
    myadata,
    basis="X_umap",
    color=color,
    frameon=False,
    ncols=1,
    save = f'_{color}_{label}.pdf'
)


output_study_paths = {"collin" : "output/scanpy/collin_merged",
"wu" : "output/scanpy/wu_merged",
"field" : "output/scanpy/field_merged",
"yang" : "output/scanpy/yang_merged"
}

metadata_paths = {
  "collin" : "data/metadata_collin_et_al.csv",
  "wu" : "data/metadata_wu_et_al.csv",
  "field" : "data/metadata_field_et_al.csv",
  "yang" : "data/metadata_yang_et_al.csv"
}

metadata = dict()
for k in metadata_paths.keys():
  # print(metadata[k])
  metadata[k] = pd.read_csv(metadata_paths[k], on_bad_lines='skip')

study_scanpy_dirs = {
  "collin" : "output/scanpy/collin_merged/",
  "yang" : "output/scanpy/yang_merged/",
  "field" : "output/scanpy/field_merged/",
  "wu" : "output/scanpy/wu_merged/"
}

metadata  = pd.read_csv("results/cell_clusters.csv")

cluster_dictionary = pd.read_csv("data/cluster_dictionary.csv")
cluster_dictionary['cluster'] = cluster_dictionary['gene_snn_res.0.2'].astype(np.float64)

cluster_dictionary = cluster_dictionary[['sample_id', 'cluster', 'abbreviation']]

metadata = metadata.merge(cluster_dictionary, on = ['sample_id', 'cluster'])

metadata.cluster = metadata['sample_id'].astype(str) + '_' + metadata['abbreviation'].astype(str)

adata_paths = {
  "SRR13884242" : "output/scanpy/SRR13884242.h5ad",
  "SRR13884243" : "output/scanpy/SRR13884243.h5ad",
  "SRR13884247" : "output/scanpy/SRR13884247.h5ad",
  "SRR13884249" : "output/scanpy/SRR13884249.h5ad",
  "SRR14800534" : "output/scanpy/SRR14800534.h5ad",
  "SRR14800535" : "output/scanpy/SRR14800535.h5ad",
  "SRR14800536" : "output/scanpy/SRR14800536.h5ad",
  "SRR14800540" : "output/scanpy/SRR14800540.h5ad",
  "SRR14800541" : "output/scanpy/SRR14800541.h5ad",
  "SRR14800543" : "output/scanpy/SRR14800543.h5ad",
  "SRR17960481" : "output/scanpy/SRR17960481.h5ad",
  "SRR17960484" : "output/scanpy/SRR17960484.h5ad"
}

def plot_dend(sample_id, adata):
  adata.obs.cluster = adata.obs.cluster.astype('category')
  sc.tl.dendrogram(adata, groupby = "cluster")
  dend = sc.pl.dendrogram(adata, groupby = "cluster", orientation = "left")
  dend.figure.set_size_inches(18, 13)
  dend.figure.savefig(f'{sample_id}_dend.pdf', dpi=200)

# need to add cluster markers from seurat objects 

adatas = dict()
for k,v in adata_paths.items():
  adatas[k] = sc.read_h5ad(adata_paths[k])
  sample_meta = metadata.loc[metadata['sample_id'] == k]
  sample_meta = sample_meta[sample_meta.cell.isin(adatas[k].obs.index)]
  sample_meta.index = sample_meta.cell
  adatas[k].obs['cluster'] = sample_meta['cluster']
  adatas[k].obs.cluster = adatas[k].obs.cluster.astype('category')

for k,v in adata_paths.items():
  # mitochondrial genes
  adatas[k].var['mt'] = adatas[k].var_names.str.startswith('MT-') 
  # ribosomal genes
  adatas[k].var['ribo'] = adatas[k].var_names.str.startswith(("RPS","RPL"))
  # hemoglobin genes.
  adatas[k].var['hb'] = adatas[k].var_names.str.contains(("^HB[^(P)]"))
  sc.pp.calculate_qc_metrics(adatas[k], qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
  sc.pp.filter_cells(adatas[k], min_genes=200)
  sc.pp.filter_genes(adatas[k], min_cells=3)
  # filter for percent mito
  adatas[k] = adatas[k][adatas[k].obs['pct_counts_mt'] < 20, :]
  # filter for percent ribo > 0.05
  adatas[k] = adatas[k][adatas[k].obs['pct_counts_ribo'] > 5, :]
  malat1 = adatas[k].var_names.str.startswith('MALAT1')
  # we need to redefine the mito_genes since they were first 
  # calculated on the full object before removing low expressed genes.
  mito_genes = adatas[k].var_names.str.startswith('MT-')
  hb_genes = adatas[k].var_names.str.contains('^HB[^(P)]')
  
  remove = np.add(mito_genes, malat1)
  remove = np.add(remove, hb_genes)
  keep = np.invert(remove)
  
  adatas[k] = adatas[k][:,keep]



# for k,v in adata_paths.items():
#   sc.pp.neighbors(adatas[k])
#   sc.pp.pca(adatas[k])
#   sc.tl.umap(adatas[k])
#   sc.tl.leiden(adatas[k], resolution = 0.2)
#   sc.tl.rank_genes_groups(adatas[k], 'leiden', method='wilcoxon')
#   sc.pl.rank_genes_groups(adatas[k], n_genes=25, sharey=False, save = f'_leiden_{k}.pdf')
#   plot_my_embedding(adatas[k], ["leiden", "cluster"], k)
#   sc.tl.dendrogram(adatas[k], groupby = "cluster")
#   dend = sc.pl.dendrogram(adatas[k], groupby = "cluster", orientation = "left")
#   dend.figure.set_size_inches(18, 13)
#   dend.figure.savefig(f'{k}_dend.pdf', dpi=200)

combined_adata = anndata.concat(adatas, label="sample_id")


combined_adata.obs_names_make_unique()

combined_adata.raw = combined_adata  # keep full dimension safe

# add numbat metadata ------------------------------

from pathlib import Path

seu_meta_paths = Path('output/seurat').glob('*embeddings*')

seu_metas = dict()
for i in seu_meta_paths:
  k = i.name.replace("_embeddings.csv", "")
  seu_metas[k] = pd.read_csv(i)
  seu_metas[k]["sample_id"] = k
  

seu_meta_df = pd.concat([v for k,v in seu_metas.items()])

seu_meta_df['cell_short'] = seu_meta_df.cell.str.replace("\\-1", "")

combined_adata.obs['cell_short'] = combined_adata.obs.index.str.replace("\\-1", "")


seu_meta_df = seu_meta_df.merge(combined_adata.obs, how = "right", on = ["cell_short", "sample_id"])

seu_meta_df.index = combined_adata.obs.index


combined_adata.obs['clone_opt'] = seu_meta_df.clone_opt.astype('category')

# # sc.pp.highly_variable_genes(
# #     combined_adata,
# #     flavor="seurat_v3",
# #     n_top_genes=2000,
# #     layer="counts",
# #     batch_key="sample_id",
# #     subset=True,
# # )
# 
# # sc.pp.highly_variable_genes(
# #     combined_adata
# # )

def scvi_integrate(myadata, adata_path):
  scvi.model.SCVI.setup_anndata(myadata, batch_key="sample_id")
  
  vae = scvi.model.SCVI(myadata, n_layers=2, n_latent=30, gene_likelihood="nb")
  
  vae.train()
  
  myadata.obsm["X_scVI"] = vae.get_latent_representation()
  
  sc.pp.neighbors(myadata, use_rep="X_scVI")
  sc.tl.leiden(myadata, resolution = 0.2)
  
  myadata.write_h5ad(adata_path)

# scvi all ------------------------------

scvi_integrate(combined_adata, "output/scanpy/scvi_adata.h5ad")

# scvi 16q ------------------------------

# kept_sample_ids = ['SRR13884242', 'SRR13884243', 'SRR13884247', 'SRR13884249', 'SRR14800534', 'SRR14800535', 'SRR14800536', 'SRR14800540', 'SRR14800541', 'SRR14800543', 'SRR17960481', 'SRR17960484']

kept_sample_ids = ['SRR14800534', 'SRR14800535', 'SRR14800536']

kept_clones = [1.0, 2.0]

adata_16q = combined_adata[combined_adata.obs.sample_id.isin(kept_sample_ids), :]

adata_16q = adata_16q[adata_16q.obs.clone_opt.isin(kept_clones), :].copy()

scvi_integrate(adata_16q, "output/scanpy/16q_adata.h5ad")

# scvi 16q plus 1q ------------------------------

kept_sample_ids = ['SRR14800534', 'SRR14800535', 'SRR14800536']

kept_clones = [2.0, 3.0]

adata_16q_plus_1q = combined_adata[combined_adata.obs.sample_id.isin(kept_sample_ids), :]

adata_16q_plus_1q = adata_16q_plus_1q[adata_16q_plus_1q.obs.clone_opt.isin(kept_clones), :].copy()

scvi_integrate(adata_16q_plus_1q, "output/scanpy/16q_plus_1q_adata.h5ad")

# scvi 1q_16q ------------------------------

kept_sample_ids = ['SRR13884242', 'SRR13884243', 'SRR13884249', 'SRR14800534', 'SRR14800535', 'SRR14800536', 'SRR14800543']

kept_clones = [1.0, 2.0, 3.0]

adata_1q_16q = combined_adata[combined_adata.obs.sample_id.isin(kept_sample_ids), :]

adata_1q_16q = adata_1q_16q[adata_1q_16q.obs.clone_opt.isin(kept_clones), :].copy()

scvi_integrate(adata_1q_16q, "output/scanpy/1q_16q_adata.h5ad")

# scvi 1q_6p ------------------------------

kept_sample_ids = ['SRR14800543', 'SRR17960484', 'SRR17960481', 'SRR17960484']

kept_clones = [1.0, 2.0, 3.0]

adata_1q_6p = combined_adata[combined_adata.obs.sample_id.isin(kept_sample_ids), :]

adata_1q_6p = adata_1q_6p[adata_1q_6p.obs.clone_opt.isin(kept_clones), :].copy()

scvi_integrate(adata_1q_6p, "output/scanpy/1q_6p_adata.h5ad")
