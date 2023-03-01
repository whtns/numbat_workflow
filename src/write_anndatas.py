
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Read in the count matrix into an AnnData object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5-based file format: .h5ad.

adata = sc.read_10x_mtx(
    mtx_dir,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

anndata_paths = [seu_path.replace("_seu.rds", ".h5ad") for seu_path in r.seu_paths]

anndata_paths = [seu_path.replace("seurat", "scanpy") for seu_path in anndata_paths]

anndata_paths = [Path(seu_path) for seu_path in anndata_paths]

sample_ids = [p.parts[1:7] + ("cellranger", p.parts[8].replace(".h5ad", ""), "outs", "filtered_feature_bc_matrix") for p in anndata_paths]

matrix_dirs = ["/".join(path_tuple) for path_tuple in sample_ids]

matrix_dirs = ["/"+"/".join(path_tuple) for path_tuple in sample_ids]

def write_my_anndata(mtx_dir):
  
  mtx_dir = Path(mtx_dir)
  
  adata = sc.read_10x_mtx(
    mtx_dir,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=False)                              # write a cache file for faster subsequent reading
  
  anndata_file = mtx_dir.parts[0:7] + ("scanpy", mtx_dir.parts[8]+".h5ad")
  
  anndata_file = "/".join(anndata_file)
  
  adata.write_h5ad(anndata_file)
  
  print(anndata_file)

[write_my_anndata(apath) for apath in matrix_dirs]
