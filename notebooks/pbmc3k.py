import matplotlib
matplotlib.use('Agg')
import scanpy as sc

sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, facecolor="white")

results_file = "pbmc3k.h5ad"

# Load data
adata = sc.datasets.pbmc3k()
adata.var_names_make_unique()

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Flag mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Filter dying cells
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Normalized.")

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes:", adata.var.highly_variable.sum())

# Save raw before filtering
adata.raw = adata

# Keep only highly variable genes
adata = adata[:, adata.var.highly_variable]
print("After HVG filter:", adata)


# Regress out technical variation
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver="arpack")
print("PCA done.")
print("PCA shape:", adata.obsm["X_pca"].shape)


# Build neighbor graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
print("Neighbors done.")

# UMAP
sc.tl.umap(adata)
print("UMAP done.")

# Leiden clustering
sc.tl.leiden(adata, resolution=0.9, random_state=0, flavor="igraph", n_iterations=2)
print("Clusters found:", adata.obs["leiden"].nunique())
print(adata.obs["leiden"].value_counts())


# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, "leiden", method="t-test")

# Print top 5 marker genes per cluster
import pandas as pd
marker_df = sc.get.rank_genes_groups_df(adata, group=None)
for cluster in sorted(adata.obs["leiden"].unique()):
    top = marker_df[marker_df["group"] == cluster].head(5)
    print(f"\nCluster {cluster}:")
    print(top[["names", "scores"]].to_string(index=False))

# Save everything
adata.write(results_file)
print("\nSaved to", results_file)


# Save UMAP plot as image
# Cell type annotations based on marker genes
cluster_names = {
    "0": "CD8+ T cells",
    "1": "CD4+ T cells",
    "2": "B cells",
    "3": "CD14+ Monocytes",
    "4": "NK cells",
    "5": "CD16+ Monocytes",
    "6": "Dendritic cells",
    "7": "Platelets",
}

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_names)

sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="right margin",
    title="PBMC 3k — Cell Types",
    save="_pbmc3k_annotated.png"
)
print("UMAP saved.")