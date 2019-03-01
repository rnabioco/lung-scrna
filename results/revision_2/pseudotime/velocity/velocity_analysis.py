#!/miniconda3/bin/python3

import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import velocyto as vcy
import pandas as pd
import loompy
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
#from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

proj_path = "/Users/kriemo/Projects/publication_repos/lung-scrna/results/revision_2/"
velo_path = os.path.join(proj_path, "pseudotime", "velocity")
vlm = vcy.VelocytoLoom(os.path.join(velo_path, "combined.loom"))
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

mdata = mdata[mdata.new_id.isin(list(vlm.ca["CellID"]))]

# reorder cell ids to match loom object
cids = pd.DataFrame({'new_id' : vlm.ca["CellID"]})
mdata = pd.merge(cids, mdata, how = 'left', on = 'new_id')
mdata = mdata.dropna()

#only keep cells found in seurat data
keep_idx = [x in list(mdata["new_id"]) for x in list(vlm.ca["CellID"])]
vlm.filter_cells(bool_array=keep_idx)

#add cluster annotations

vlm.ca["Cluster"] = np.array([int(x) for x in mdata["cluster"]])

#add tSNE projections
vlm.ca["TSNE1"] = np.array(mdata["tSNE_1"])
vlm.ca["TSNE2"] = np.array(mdata["tSNE_2"])
vlm.ts = np.column_stack([vlm.ca["TSNE1"], vlm.ca["TSNE2"]])


vlm.normalize("S", size=True, log=True)

vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

# get from seurat obj
vlm.set_clusters(vlm.ca["Cluster"])

vlm.score_detection_levels(min_expr_counts=30, min_cells_express=20)
vlm.filter_genes(by_detection_levels=True)
vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)

vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k = int(len(vlm.ca["CellID"]) * 0.025), n_jobs=6 )
#vlm.knn_imputation(n_pca_dims=20, k = len(vlm.ca["CellID"]) * 0.025, balanced=True, b_sight=3000,b_maxl=1500, n_jobs=6 )
vlm.fit_gammas()
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)


vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             knn_random=True, sampled_fraction=0.3,
                             random_seed = 42)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
vlm.calculate_grid_arrows(smooth=0.5, steps=(40, 40), n_neighbors=50)


vlm.to_hdf5("combined.hdf5")

###


###
vlm = vcy.load_velocyto_hdf5("combined.hdf5")


def despline():
    ax1 = plt.gca()
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

def minimal_xticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    xlims = np.linspace(start, end_, 5)
    xlims_tx = [""]*len(xlims)
    xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
    plt.xticks(xlims, xlims_tx)


def minimal_yticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    ylims = np.linspace(start, end_, 5)
    ylims_tx = [""]*len(ylims)
    ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
    plt.yticks(ylims, ylims_tx)



plt.figure(None,(14,14))
quiver_scale = 10
plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
            c="0.8", alpha=0.2, s=10, edgecolor="")

ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/1.), replace=False)
plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
            c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,linewidths=0.25,
                   width=0.00045,
                   edgecolors="k",
                   color=vlm.colorandum[ix_choice],
                   alpha=1)
plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
           vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
           scale=quiver_scale, **quiver_kwargs)

plt.axis("off")

plt.savefig("full_arrows.pdf")

plt.figure(None,(20,10))
vlm.plot_grid_arrows(quiver_scale=3.0,
                     plot_random=True,
                     scale_type="relative")
plt.savefig("vectorfield.pdf")

genes = ["Pdpn","Hopx", "Emp2", "Trp53", "Top2a", "Aqp5", "Rtkn2", "Ager"]
plt.figure(None, (17,24), dpi=80)
gs = plt.GridSpec(10,6)
for i, gn in enumerate(genes):
    ax = plt.subplot(gs[i*3])
    try:
        ix=np.where(vlm.ra["Gene"] == gn)[0][0]
    except:
        continue
    vcy.scatter_viz(vlm.Sx_sz[ix,:], vlm.Ux_sz[ix,:], c=vlm.colorandum, s=5, alpha=0.4, rasterized=True)
    plt.title(gn)
    xnew = np.linspace(0,vlm.Sx[ix,:].max())
    plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
    plt.ylim(0, np.max(vlm.Ux_sz[ix,:])*1.02)
    plt.xlim(0, np.max(vlm.Sx_sz[ix,:])*1.02)
    minimal_yticks(0, np.max(vlm.Ux_sz[ix,:])*1.02)
    minimal_xticks(0, np.max(vlm.Sx_sz[ix,:])*1.02)
    despline()

    vlm.plot_velocity_as_color(gene_name=gn, gs=gs[i*3+1], s=3, rasterized=True)

    vlm.plot_expression_as_color(gene_name=gn, gs=gs[i*3+2], s=3, rasterized=True)

plt.tight_layout()
plt.savefig("genes.pdf")


## Markov Chain analysis

# Sample uniformly the points to avoid density driven effects - Should reimplement as a method
steps = 100, 100
grs = []
for dim_i in range(vlm.embedding.shape[1]):
    m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
    m = m - 0.025 * np.abs(M - m)
    M = M + 0.025 * np.abs(M - m)
    gr = np.linspace(m, M, steps[dim_i])
    grs.append(gr)

meshes_tuple = np.meshgrid(*grs)
gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

from sklearn.neighbors import NearestNeighbors
nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(gridpoints_coordinates, 1)

diag_step_dist = np.sqrt((meshes_tuple[0][0,0] - meshes_tuple[0][0,1])**2 + (meshes_tuple[1][0,0] - meshes_tuple[1][1,0])**2)
min_dist = diag_step_dist / 2
ixs = ixs[dist < min_dist]
gridpoints_coordinates = gridpoints_coordinates[dist.flat[:]<min_dist,:]
dist = dist[dist < min_dist]

ixs = np.unique(ixs)

# plt.figure(None,(8,8))
# vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
#                 c=vlm.colorandum[ixs], alpha=1, s=30, lw=0.4,
#                 edgecolor="0.4")
# plt.savefig("clusters.pdf")

vlm.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2., direction='forward', cells_ixs=ixs)

vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=5000)


diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
diffused_n /= np.percentile(diffused_n, 97)
diffused_n = np.clip(diffused_n, 0, 1)

plt.figure(None,(7,7))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.5, s=50, lw=0.,
                edgecolor="", cmap="viridis_r", rasterized=False)
plt.axis("off")
plt.savefig("endpoint.pdf")


vlm.prepare_markov(sigma_D=diag_step_dist, sigma_W=diag_step_dist/2., direction='backwards', cells_ixs=ixs)
vlm.run_markov(starting_p=np.ones(len(ixs)), n_steps=5000)

diffused_n = vlm.diffused - np.percentile(vlm.diffused, 3)
diffused_n /= np.percentile(diffused_n, 97)
diffused_n = np.clip(diffused_n, 0, 1)

plt.figure(None,(7,7))
vcy.scatter_viz(vlm.embedding[ixs, 0], vlm.embedding[ixs, 1],
                c=diffused_n, alpha=0.5, s=50, lw=0.,
                edgecolor="", cmap="viridis_r", rasterized=False)
plt.axis("off")
plt.savefig("startpoint.pdf")






## examine transistion probabilities
def gaussian_kernel(X, mu = 0, sigma=1):
    return np.exp(-(X - mu)**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)

from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors

plt.figure(None,(6,6))

steps = 45, 45
grs = []
for dim_i in range(vlm.embedding.shape[1]):
    m, M = np.min(vlm.embedding[:, dim_i]), np.max(vlm.embedding[:, dim_i])
    m = m - 0.025 * np.abs(M - m)
    M = M + 0.025 * np.abs(M - m)
    gr = np.linspace(m, M, steps[dim_i])
    grs.append(gr)

meshes_tuple = np.meshgrid(*grs)
gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T
gridpoints_coordinates = gridpoints_coordinates + norm.rvs(loc=0, scale=0.15, size=gridpoints_coordinates.shape)

nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(gridpoints_coordinates, 50)
ix_choice = ixs[:,0].flat[:]
ix_choice = np.unique(ix_choice)

nn = NearestNeighbors()
nn.fit(vlm.embedding)
dist, ixs = nn.kneighbors(vlm.embedding[ix_choice], 50)
density_extimate = gaussian_kernel(dist, mu=0, sigma=0.5).sum(1)
bool_density = density_extimate > np.percentile(density_extimate, 25)
ix_choice = ix_choice[bool_density]

plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
            c=vlm.colorandum, alpha=0.2, s=120, edgecolor="")
plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
            c=vlm.colorandum[ix_choice], alpha=1, s=120, edgecolor="k")

quiver_kwargs=dict(scale=1, headaxislength=9, headlength=15, headwidth=14,linewidths=0.4, edgecolors="k", color="k", alpha=1)
plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
           vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
           **quiver_kwargs)

plt.xlim(-50.,-20)
plt.ylim(-20,10)
plt.savefig("grid_arrows_cca.pdf")
