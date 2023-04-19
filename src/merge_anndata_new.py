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


def scvi_combined_adata(combined_adata, processed_study_path):
  
  sc.pp.highly_variable_genes(
      combined_adata,
      flavor="seurat_v3",
      n_top_genes=2000,
      batch_key="sample_id",
      subset=True
  )
  
  sc.tl.pca(combined_adata)
  
  scvi.model.SCVI.setup_anndata(combined_adata, batch_key="sample_id")

  vae = scvi.model.SCVI(combined_adata, n_layers=2, n_latent=30, gene_likelihood="nb")

  vae.train()

  combined_adata.obsm["X_scVI"] = vae.get_latent_representation()

  sc.pp.neighbors(combined_adata, use_rep="X_scVI")

  combined_adata.obsm["X_mde"] = mde(combined_adata.obsm["X_scVI"])
  
  combined_adata.write_h5ad(processed_study_path)

  # sc.pl.embedding(
  #     cross_study_adata,
  #     basis="X_mde",
  #     color=["batch", "leiden"],
  #     frameon=False,
  #     ncols=1,
  # )
  
  # sc.external.pp.bbknn(combined_adata, batch_key='sample_id') 
  
  # sc.pp.pca(combined_adata)
  # sc.pp.neighbors(combined_adata)
  sc.tl.umap(combined_adata)
  sc.tl.leiden(combined_adata)
  
  
  return combined_adata


def integrate_combined_adata(combined_adata, processed_study_path):
  
  sc.pp.highly_variable_genes(
      combined_adata,
      flavor="seurat_v3",
      n_top_genes=2000,
      batch_key="sample_id",
      subset=True
  )
  
  sc.tl.pca(combined_adata)
  
  sc.external.pp.bbknn(combined_adata, batch_key='sample_id') 
  
  # sc.pp.pca(combined_adata)
  # sc.pp.neighbors(combined_adata)
  sc.tl.umap(combined_adata)
  sc.tl.leiden(combined_adata)
  
  combined_adata.write_h5ad(processed_study_path)
  
  return combined_adata


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
  "SRR13884242" : "output/scanpy/wu_merged/SRR13884242.h5ad",
  "SRR13884243" : "output/scanpy/wu_merged/SRR13884243.h5ad",
  "SRR13884247" : "output/scanpy/wu_merged/SRR13884247.h5ad",
  "SRR13884249" : "output/scanpy/wu_merged/SRR13884249.h5ad",
  "SRR14800534" : "output/scanpy/yang_merged/SRR14800534.h5ad",
  "SRR14800535" : "output/scanpy/yang_merged/SRR14800535.h5ad",
  "SRR14800536" : "output/scanpy/yang_merged/SRR14800536.h5ad",
  "SRR14800540" : "output/scanpy/yang_merged/SRR14800540.h5ad",
  "SRR14800541" : "output/scanpy/yang_merged/SRR14800541.h5ad",
  "SRR14800543" : "output/scanpy/yang_merged/SRR14800543.h5ad",
  "SRR17960481" : "output/scanpy/field_merged/SRR17960481.h5ad",
  "SRR17960484" : "output/scanpy/field_merged/SRR17960484.h5ad"
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
  sc.tl.dendrogram(adatas[k], groupby = "cluster")
  dend = sc.pl.dendrogram(adatas[k], groupby = "cluster", orientation = "left")
  dend.figure.set_size_inches(18, 13)
  dend.figure.savefig(f'{k}_dend.pdf', dpi=200)


combined_adata = anndata.concat(adatas, label="sample_id")

# combined_adata0 = combined_adata.copy()
# combined_adata0 = integrate_combined_adata(combined_adata0, "output/scanpy/interesting_adata.h5ad")
combined_adata0 = sc.read_h5ad("output/scanpy/interesting_adata.h5ad")

# # combined_adata0.obs.cluster = combined_adata0.obs.cluster.astype(str)

# sc.pl.umap(combined_adata0, color = ['leiden', 'sample_id', 'cluster'])
# sc.tl.leiden(combined_adata0, key_added = "leiden_1.0")

sc.tl.dendrogram(combined_adata0, groupby = "cluster")
dend_integrated = sc.pl.dendrogram(combined_adata0, groupby = "cluster", orientation = "left")

dend_integrated.figure.set_size_inches(18, 13)
dend_integrated.figure.savefig("integrated_umap.pdf", dpi=200)

# sc.tl.pca(combined_adata)
combined_adata.obs.cluster = combined_adata.obs.cluster.astype('category')
sc.tl.dendrogram(combined_adata, groupby = "cluster")
dend = sc.pl.dendrogram(combined_adata, groupby = "cluster", orientation = "left")

dend.figure.set_size_inches(18, 13)
dend.figure.savefig("umap.pdf", dpi=200)


sc.pp.log1p(combined_adata)

sc.tl.rank_genes_groups(combined_adata, 'cluster', method='wilcoxon')

# sc.pl.rank_genes_groups_heatmap(combined_adata, save = "groups_heatmap.pdf", key = "cluster", show_gene_labels=True)

# first gene
sc.pl.rank_genes_groups_matrixplot(
    combined_adata,
    n_genes=1,
    values_to_plot="logfoldchanges",
    cmap='bwr',
    vmin=-4,
    vmax=4,
    min_logfoldchange=3,
    colorbar_title='log fold change',
    save = "matrix_plot"
)

# last gene
sc.pl.rank_genes_groups_matrixplot(
    combined_adata,
    n_genes=-1,
    values_to_plot="logfoldchanges",
    cmap='bwr',
    vmin=-4,
    vmax=4,
    min_logfoldchange=3,
    colorbar_title='log fold change',
    save = "matrix_plot"
)


markers = ['MARCKSL1', 'MEG3', 'ARR3', 'TOP2A', 'B2M']
matrixplot = sc.pl.matrixplot(combined_adata, markers, groupby = "cluster", save = "matrixplot.pdf", dendrogram=True)

matrixplot.figure.set_size_inches(18, 13)
matrixplot.figure.savefig("matrixplot.pdf", dpi=200)




# process merged study adatas by bbknn ------------------------------

def process_study_adata(study_path):

  study_adata = sc.read_h5ad(study_path)
  
  processed_study_path = str(study_path).replace(".h5ad", "_processed.h5ad")
  
  sc.pp.highly_variable_genes(
      study_adata,
      flavor="seurat_v3",
      n_top_genes=2000,
      batch_key="sample_id",
      subset=True
  )
  
  sc.tl.pca(study_adata)
  
  sc.external.pp.bbknn(study_adata, batch_key='sample_id') 
  
  # sc.pp.pca(study_adata)
  # sc.pp.neighbors(study_adata)
  sc.tl.umap(study_adata)
  sc.tl.leiden(study_adata)
  
  print(processed_study_path)
  study_adata.write_h5ad(processed_study_path)

for study in merged_study_paths:
  process_study_adata(merged_study_paths[study])

# collin ------------------------------

collin_processed_adata = sc.read_h5ad("/home/skevin/single_cell_projects/resources/collin_et_al_proj/output/scanpy/merged/merged_study_processed.h5ad")

sc.pl.umap(processed_collin_adata, color = ['leiden', 'sample_id'])

sc.tl.leiden(collin_processed_adata, resolution = 0.2)
sc.pl.umap(collin_processed_adata, color = ['leiden', 'sample_id', 'TOP2A'], save = ".png")

sc.pp.log1p(collin_processed_adata)

sc.tl.leiden(collin_processed_adata, resolution = 0.2)
sc.pl.umap(collin_processed_adata, color = ['leiden', 'sample_id'], save = ".png")

sc.tl.rank_genes_groups(collin_processed_adata, 'leiden', method='wilcoxon')

# to visualize the results
sc.pl.rank_genes_groups(collin_processed_adata, n_genes = 10)

sc.pl.umap(collin_processed_adata, color = ['leiden', 'MAP1B', 'KRT13', 'VIM', 'HIST1H4C', 'HBA2'], save = ".png")

# field ------------------------------

field_processed_adata = sc.read_h5ad("/home/skevin/single_cell_projects/resources/field_et_al_proj/output/scanpy/merged/merged_study_processed.h5ad")

sc.pp.log1p(field_processed_adata)

sc.tl.leiden(field_processed_adata, resolution = 0.2)
sc.pl.umap(field_processed_adata, color = ['leiden', 'sample_id'], save = ".png")

sc.tl.rank_genes_groups(field_processed_adata, 'leiden', method='wilcoxon')

# to visualize the results
sc.pl.rank_genes_groups(field_processed_adata, n_genes = 10)

sc.pl.umap(field_processed_adata, color = ['leiden', 'GABARAP', 'IGFBP3', 'HIST1H4C', 'S100A11'], save = ".png")

study = "field"
merged_metada_path = f"~/single_cell_projects/resources/{study}_et_al_proj/output/scanpy/merged/metadata.csv"

field_processed_adata.obs.to_csv(merged_metada_path)

# scvi ------------------------------

# scvi.model.SCVI.setup_anndata(collin_adata, batch_key="sample_id")
# 
# vae = scvi.model.SCVI(collin_adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# 
# vae.train()
# 
# cross_study_adata.obsm["X_scVI"] = vae.get_latent_representation()
# 
# sc.pp.neighbors(cross_study_adata, use_rep="X_scVI")
# sc.tl.leiden(cross_study_adata)
# 
# from scvi.model.utils import mde
# import pymde
# cross_study_adata.obsm["X_mde"] = mde(cross_study_adata.obsm["X_scVI"])
# 
# sc.pl.embedding(
#     cross_study_adata,
#     basis="X_mde",
#     color=["batch", "leiden"],
#     frameon=False,
#     ncols=1,
# )

# yang ------------------------------

yang_processed_adata = sc.read_h5ad("/home/skevin/single_cell_projects/resources/yang_et_al_proj/output/scanpy/merged/merged_study_processed.h5ad")

sc.pl.umap(yang_processed_adata, color = ['leiden', 'sample_id', 'TOP2A'], save = ".png")

marker_genes_dict = {
  'cone': ['NRL'],
  'lymph': 'B2M'}

ax = sc.pl.dotplot(yang_processed_adata, marker_genes_dict, groupby='leiden', dendrogram=True,
                   standard_scale='var', smallest_dot=40, color_map='Blues', figsize=(8,5))

gs = sc.pl.matrixplot(yang_processed_adata, marker_genes_dict, groupby='leiden')

sc.pl.umap(yang_processed_adata, color = ['B2M', 'NRL'], save = ".png")

sc.pp.log1p(yang_processed_adata)

sc.tl.leiden(yang_processed_adata, resolution = 0.2)
sc.pl.umap(yang_processed_adata, color = ['leiden', 'sample_id'], save = ".png")

sc.tl.rank_genes_groups(yang_processed_adata, 'leiden', method='wilcoxon')

# to visualize the results
sc.pl.rank_genes_groups(yang_processed_adata, n_genes = 10)

sc.pl.umap(yang_processed_adata, color = ['leiden', 'TFF1', 'MT-CO1', 'LGALS3', 'VIM', 'TOP2A', 'ROM1', 'CD52'], save = ".png")

yang_processed_adata.obs.to_csv("/home/skevin/single_cell_projects/resources/yang_et_al_proj/output/scanpy/merged/metadata.csv")


# wu ------------------------------

wu_processed_adata = sc.read_h5ad("output/scanpy/wu_merged/merged_study_processed.h5ad")

sc.pp.log1p(wu_processed_adata)

sc.tl.leiden(wu_processed_adata, resolution = 0.2)
sc.pl.umap(wu_processed_adata, color = ['leiden', 'sample_id'], save = ".png")

sc.tl.rank_genes_groups(wu_processed_adata, 'leiden', method='wilcoxon')

# to visualize the results

sc.pl.rank_genes_groups(wu_processed_adata, n_genes = 10)

sc.pl.umap(wu_processed_adata, color = ['leiden', 'RPL6', 'GNAT1', 'HIST1H4C', 'ARL6IP1', 'FTL', 'GUCA1B'], save = ".png")

wu_processed_adata.obs.to_csv("/home/skevin/single_cell_projects/resources/wu_et_al_proj/output/scanpy/merged/metadata.csv")


# cross study ------------------------------
cross_study_adata = anndata.concat(study_adatas, label = "study")

cross_study_adata.raw = cross_study_adata  # keep full dimension safe

# cross_study_adata.write_h5ad("output/scanpy/cross_study_adata.h5ad")
cross_study_adata = sc.read_h5ad("output/scanpy/cross_study_adata.h5ad")

subset_cross_study_adata = cross_study_adata

metadata = pd.read_csv("results/cell_clusters.csv")



integrate_combined_adata("output/scanpy/cross_study_adata.h5ad")

integrated_combined_adata = sc.read_h5ad("output/scanpy/cross_study_adata_bbknn_integrated.h5ad")

sc.pl.umap(integrated_combined_adata, color = ['leiden', 'sample_id', 'study', 'TOP2A', 'NRL', 'RCVRN'], save = ".png")

sc.pl.umap(integrated_combined_adata, color = ['sample_id'], legend_fontsize = "xx-large", save = ".png")


# scvi ------------------------------
scvi.model.SCVI.setup_anndata(cross_study_adata, batch_key="sample_id")

vae = scvi.model.SCVI(cross_study_adata, n_layers=2, n_latent=30, gene_likelihood="nb")

vae.train()

cross_study_adata.obsm["X_scVI"] = vae.get_latent_representation()

sc.pp.neighbors(cross_study_adata, use_rep="X_scVI")
sc.tl.leiden(cross_study_adata)

from scvi.model.utils import mde
import pymde
cross_study_adata.obsm["X_mde"] = mde(cross_study_adata.obsm["X_scVI"])

sc.pl.embedding(
    cross_study_adata,
    basis="X_mde",
    color=["batch", "leiden"],
    frameon=False,
    ncols=1,
)
# ------------------------------

sc.pp.pca(cross_study_adata)
sc.pp.neighbors(cross_study_adata)
sc.tl.umap(cross_study_adata)

sc.pl.umap(cross_study_adata, color=['study', 'TFF1', 'NRL', 'HIST1H4C'], legend_loc = "upper center")


cross_study_adata.write_h5ad("output/scanpy/cross_study_adata_processed.h5ad")

sc.tl.pca(cross_study_adata, svd_solver='arpack')

sc.pl.pca(cross_study_adata, dimensions = [(5,6)], color='RXRG')

# read adata ------------------------------
cross_study_adata = sc.read_h5ad("output/scanpy/subset_adata.h5ad")


cross_study_adata.obs['combined'] = cross_study_adata.obs['cell'].astype(str) + '_' + cross_study_adata.obs['sample_id'].astype(str)
cross_study_adata.obs.index = cross_study_adata.obs['combined']

metadata = pd.read_csv("results/cell_clusters.csv")

metadata['combined'] = metadata['cell'].astype(str) + '_' + metadata['sample_id'].astype(str)
metadata.index = metadata.combined

metadata = metadata[metadata.combined.isin(cross_study_adata.obs.combined)]

cross_study_adata2 = cross_study_adata.copy()

cross_study_adata0 = cross_study_adata2[cross_study_adata2.obs.combined.isin(metadata.combined)]


# test0 ------------------------------





