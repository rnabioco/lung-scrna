#!/miniconda3/bin/python3

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
velo_path = os.path.join(proj_path, "pseudotime", "velocity")

adata = scv.read_loom(os.path.join(velo_path, "combined.loom"), sparse = True)
adata.var_names_make_unique()

mdata = pd.read_csv(os.path.join(proj_path, "revision/geo/cell_metadata.tsv.gz"), sep="\t")

sample_map = pd.read_csv(os.path.join(velo_path, "sample_to_path.tsv"), sep = "\t")

sample = [ x.split("/")[8] for x in sample_map['paths'].tolist()]
new_id = sample_map['samples'].tolist()

good_clusters = [
  'Naive AEC1',
  "Injured AEC2: Cell Cycle Arrest",
  "Injured AEC2: Proliferating",
  "Other Injured AEC2",
  "Naive AEC2",
  "Injured AEC2: Transdifferentiating"]


mdata = mdata[mdata["cell_type"].isin(good_clusters)]

sample_id_map = {}
for i in range(len(new_id)):
    sample_id_map[new_id[i]] = sample[i]

# get cell ids to match loom object
new_ids  = [x.replace(x.split("_")[0] + "_", sample_id_map[x.split("_")[0]] + ":") + "x" for x in mdata["cell"] ]

mdata = mdata.assign(new_id = new_ids)

mdata = mdata[mdata.new_id.isin(list(adata.obs.index.values))]

# reorder cell ids to match loom object
cids = pd.DataFrame({'new_id' : adata.obs.index.values})
mdata = pd.merge(cids, mdata, how = 'left', on = 'new_id')
mdata = mdata.dropna()

#only keep cells found in seurat data
keep_idx = [x in list(mdata["new_id"]) for x in list(adata.obs.index.values)]

adata = adata[keep_idx, :]
adata = adata.copy()

#add cluster annotations

adata.obs['clusters'] = np.array([str(x) for x in mdata["cluster"]])
adata.obs['cell_type'] = np.array(mdata["cell_type"])

#add tSNE projections
tsne_mat = np.column_stack((np.array(mdata["tSNE_1"]), np.array(mdata["tSNE_2"])))
adata.obsm['X_tsne'] = tsne_mat

scv.utils.show_proportions(adata)

scv.pp.filter_genes(adata, min_counts=20, min_counts_u=10)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=3000)
scv.pp.log1p(adata)


scv.pp.moments(adata, n_pcs=20, n_neighbors=30)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.tl.velocity_embedding(adata, basis='tsne')
scv.pl.velocity_embedding_grid(adata)
scv.pl.velocity_embedding_stream(adata, legend_loc='on data', color = 'cell_type')
scv.pl.velocity_embedding(adata, scale = 10, basis='tsne', dpi=200)

scv.tl.cell_fate(adata)
scv.tl.terminal_states(adata)
scv.pl.velocity_embedding_grid(adata, scale = 0.1, color=['cell_fate', 'root_cells', 'end_points'], legend_loc='on data', dpi = 300)

scv.tl.velocity_confidence(adata)

scv.tl.paga(adata, groups='clusters')
scv.tl.paga(adata, groups='clusters', use_rna_velocity=True)
scv.pl.paga_compare(adata, basis='umap', threshold=.15, arrowsize=10, edge_width_scale=.5,
                    transitions='transitions_confidence', dashed_edges='connectivities')
