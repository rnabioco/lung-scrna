import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc
import pandas as pd
import loompy

scv.settings.set_figure_params('scvelo')

proj_path = "/Users/kriemo/Projects/publication_repos/lung-scrna/results/revision_2/"
input_path = os.path.join(proj_path, "revision", "geo")

adata = sc.read_text(os.path.join(input_path, "count_matrix.tsv.gz"), delimiter= "\t")
adata = adata.T

mdata = pd.read_csv(os.path.join(proj_path, "revision/geo/cell_metadata.tsv.gz"), sep="\t")

adata.obs["cluster"] = np.array(mdata["cluster"])
adata.obs["cell_type"] = np.array(mdata["cell_type"])

#add tSNE projections
tsne_mat = np.column_stack((np.array(mdata["tSNE_1"]), np.array(mdata["tSNE_2"])))
adata.obsm['X_tsne'] = tsne_mat

good_clusters = [
  'Naive AEC1',
  "Injured AEC2: Cell Cycle Arrest",
  "Injured AEC2: Proliferating",
  "Other Injured AEC2",
  "Naive AEC2",
  "Injured AEC2: Transdifferentiating"]

adata = adata[np.isin(adata.obs["cell_type"].values, np.array(good_clusters)), :]

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)


sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata, layout = "fa")
sc.pl.draw_graph(adata, color='cell_type', legend_loc='on data')

sc.tl.paga(adata, groups='cell_type')

sc.tl.draw_graph(adata, init_pos='paga')

sc.pl.paga_compare(
    adata, threshold=0.01, title='', right_margin=0, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)
