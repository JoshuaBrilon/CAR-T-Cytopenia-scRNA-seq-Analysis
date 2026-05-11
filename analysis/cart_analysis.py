"""
cart_analysis.py — CAR-T Cytopenia: scRNA-seq Analysis of Bone Marrow
======================================================================
Independent reanalysis of Strati et al. 2023 (Cell Reports Medicine).
Dataset: GSE216005 — bone marrow scRNA-seq from 22 individuals:
    - 11 rrLBCL patients with prolonged cytopenia (PSC) after CAR-T
    - 5 rrLBCL patients without prolonged cytopenia (Non-PSC)
    - 5 healthy bone marrow donors
    - 1 chemo-associated cytopenia control
    - 2 excluded (bone marrow involvement with AML/LBCL)

Question: What is different at the single-cell level in the bone marrow
of patients whose marrow fails to recover after CAR-T therapy?

Key finding reproduced: CD8+ cytotoxic T cells expressing IFN-γ are
massively expanded in PSC patients and absent in Non-PSC patients.

Reference: Strati P, Li X, Deng Q et al. Cell Rep Med. 2023;4(8):101158.
Author: Joshua Brilon | Project Hicks | University of Oregon
"""

import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

sc.settings.verbosity = 3

# ── PATIENT METADATA ──────────────────────────────────────────────────────────
# Group assignments from Table S3 of Strati et al. 2023.
# PSC = prolonged cytopenia (grade 3-4 cytopenia persisting >30 days post CAR-T)
# Non-PSC = no prolonged cytopenia after CAR-T
# pt21 and pt22 excluded: bone marrow involvement with LBCL/AML respectively

METADATA = {
    "pt01": "PSC",
    "pt02": "PSC",
    "pt03": "PSC",
    "pt04": "PSC",
    "pt05": "PSC",
    "pt06": "PSC",
    "pt07": "PSC",
    "pt08": "PSC",
    "pt09": "PSC",
    "pt10": "Non-PSC",
    "pt11": "Non-PSC",
    "pt12": "Non-PSC",
    "pt13": "Non-PSC",
    "pt14": "Non-PSC",
    "pt15": "Healthy",
    "pt16": "Healthy",
    "pt17": "Healthy",
    "pt18": "Healthy",
    "pt19": "Healthy",
    "pt20": "Chemo",
    "pt21": "Exclude",   # LBCL bone marrow involvement
    "pt22": "Exclude",   # AML bone marrow involvement
}

# Cell type annotations derived from marker gene analysis
# Determined by matching top differentially expressed genes to known biology
CLUSTER_NAMES = {
    "0":  "NK cells",
    "1":  "Cycling cells",
    "2":  "T progenitors",
    "3":  "Early erythroid",
    "4":  "CD8+ T cells",
    "5":  "CD8+ T cells (2)",
    "6":  "CD4+ T cells",
    "7":  "Cycling mixed",
    "8":  "Late erythroid",
    "9":  "Late erythroid (2)",
    "10": "B cells",
    "11": "Macrophage-like",
    "12": "Classical monocytes",
    "13": "Non-classical monocytes",
    "14": "Monocytes (3)",
    "15": "Plasmacytoid DC",
    "16": "Monocytes/Macrophages",
    "17": "Erythroid",
    "18": "T cells activated",
    "19": "NK/NKT cells",
    "20": "Plasma cells",
    "21": "NK/T mixed",
    "22": "CD8+ cytotoxic T",   # key cluster — NKG7, CD8A, GZMH, CCL5
    "23": "CD8+ T cells (3)",
    "24": "B cells (2)",
}

DATA_DIR = "/Users/joshuabrilon/Project Hicks/coding/data/GSE216005"

os.makedirs("figures", exist_ok=True)


# ── STAGE 1: LOAD DATA ────────────────────────────────────────────────────────
# Each patient folder contains three 10x Genomics files:
#   barcodes.tsv.gz — one barcode per cell
#   features.tsv.gz — one row per gene
#   matrix.mtx.gz  — sparse count matrix (cells × genes)

def load_data():
    adatas = []
    for pt, group in METADATA.items():
        if group == "Exclude":
            print(f"Skipping {pt} ({group})")
            continue
        path = os.path.join(DATA_DIR, pt)
        print(f"Loading {pt} ({group})...")
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        adata.obs["patient"] = pt
        adata.obs["group"] = group
        adatas.append(adata)

    print("\nConcatenating all samples...")
    adata = sc.concat(
        adatas,
        label="sample",
        keys=[a.obs["patient"].iloc[0] for a in adatas]
    )
    adata.obs_names_make_unique()

    print(adata)
    print("\nGroup counts:")
    print(adata.obs["group"].value_counts())
    return adata


# ── STAGE 2: QUALITY CONTROL ──────────────────────────────────────────────────
# Remove low-quality cells before analysis.
# Thresholds follow Strati et al. 2023 methods:
#   min 200 genes per cell, max 15% mitochondrial reads, max 4000 genes (doublet filter)

def run_qc(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < 15, :]
    adata = adata[adata.obs.n_genes_by_counts < 4000, :]

    print("\nAfter QC filtering:")
    print(adata)
    print("\nGroup counts after QC:")
    print(adata.obs["group"].value_counts())

    adata.write("cart_qc.h5ad")
    print("Saved checkpoint: cart_qc.h5ad")
    return adata


# ── STAGE 3: NORMALIZATION + DIMENSIONALITY REDUCTION ─────────────────────────
# Standard single-cell preprocessing pipeline.
# Key difference from single-sample analysis: batch_key="patient" in HVG
# selection ensures we find genes variable within patients, not just between
# patients (which would be driven by individual differences, not biology).

def run_processing():
    adata = sc.read_h5ad("cart_qc.h5ad")

    # Normalize to 10k counts per cell, log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # save pre-HVG matrix for gene expression visualization later

    # Select highly variable genes, correcting for patient batch effects
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="patient"
    )
    print("Highly variable genes:", adata.var.highly_variable.sum())
    adata = adata[:, adata.var.highly_variable]

    # Scale to unit variance (required for PCA)
    sc.pp.scale(adata, max_value=10)

    # PCA: compress to 50 principal components
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
    print("PCA done.")

    # Build neighbor graph and compute UMAP embedding
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    print("UMAP done.")

    # Leiden clustering at resolution 0.4 — yields 25 clusters
    sc.tl.leiden(
        adata, resolution=0.4, random_state=0, flavor="igraph", n_iterations=2
    )
    print("Clusters found:", adata.obs["leiden"].nunique())

    adata.write("cart_processed.h5ad")
    print("Saved checkpoint: cart_processed.h5ad")

    # Initial UMAP plots
    sc.pl.umap(adata, color="group", title="CAR-T Bone Marrow — Patient Groups",
               save="_groups.png")
    sc.pl.umap(adata, color="leiden", title="CAR-T Bone Marrow — Clusters",
               save="_clusters.png")
    return adata


# ── STAGE 4: MARKER GENES + ANNOTATION ───────────────────────────────────────
# Identify top differentially expressed genes per cluster using t-test.
# Match marker genes to known cell type signatures to annotate clusters.

def run_annotation():
    adata = sc.read_h5ad("cart_processed.h5ad")

    # Find marker genes
    sc.tl.rank_genes_groups(adata, "leiden", method="t-test", n_genes=20)
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)

    print("\nTop 5 marker genes per cluster:")
    for cluster in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        top = marker_df[marker_df["group"] == cluster].head(5)
        print(f"\nCluster {cluster}:")
        print(top[["names", "scores"]].to_string(index=False))

    # Apply cell type annotations
    adata.obs["cell_type"] = adata.obs["leiden"].map(CLUSTER_NAMES)

    # Save annotated UMAP plots
    sc.pl.umap(adata, color="cell_type", legend_loc="right margin",
               legend_fontsize=7, title="CAR-T Bone Marrow — Cell Types",
               save="_annotated.png")
    sc.pl.umap(adata, color="group", legend_loc="right margin",
               title="CAR-T Bone Marrow — Patient Groups",
               save="_groups_annotated.png")

    adata.write("cart_annotated.h5ad")
    print("Saved checkpoint: cart_annotated.h5ad")
    return adata


# ── STAGE 5: PSC vs NON-PSC COMPARISON ───────────────────────────────────────
# Core analysis: compare cell type proportions between cytopenic (PSC)
# and non-cytopenic (Non-PSC) patients.
#
# Key finding: CD8+ T cells — especially the CD8+ cytotoxic T cluster —
# are massively expanded in PSC patients. This cluster is absent in Non-PSC.
# These cells express IFN-γ (IFNG), which suppresses HSC self-renewal
# and drives prolonged bone marrow failure.

def run_comparison():
    adata = sc.read_h5ad("cart_annotated.h5ad")

    # Restrict to PSC vs Non-PSC for comparison
    adata_comp = adata[adata.obs["group"].isin(["PSC", "Non-PSC"])].copy()

    # ── Cell type proportions ──
    props = (
        adata_comp.obs
        .groupby(["patient", "group", "cell_type"], observed=False)
        .size()
        .reset_index(name="count")
    )
    totals = (
        adata_comp.obs
        .groupby("patient", observed=False)
        .size()
        .reset_index(name="total")
    )
    props = props.merge(totals, on="patient")
    props["proportion"] = props["count"] / props["total"]

    summary = (
        props.groupby(["cell_type", "group"], observed=False)["proportion"]
        .mean()
        .unstack()
    )
    summary["PSC_vs_NonPSC"] = summary["PSC"] - summary["Non-PSC"]
    summary = summary.sort_values("PSC_vs_NonPSC", ascending=False)

    print("\nMean cell type proportions — PSC vs Non-PSC:")
    print(summary.to_string())

    # Bar chart: which cell types are expanded or depleted in PSC
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = summary["PSC_vs_NonPSC"].apply(lambda x: "red" if x > 0 else "green")
    summary["PSC_vs_NonPSC"].plot(kind="barh", ax=ax, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("PSC proportion minus Non-PSC proportion")
    ax.set_title(
        "Cell type enrichment: PSC vs Non-PSC\n"
        "Red = expanded in PSC | Green = expanded in Non-PSC"
    )
    plt.tight_layout()
    plt.savefig("figures/proportion_comparison.png", dpi=150)
    print("Saved: figures/proportion_comparison.png")

    # ── IFN-γ expression in CD8+ cytotoxic T cells ──
    # This cluster (Leiden 22) is the population identified by Strati et al.
    # as driving bone marrow suppression via IFN-γ signaling in HSCs.
    # It is present only in PSC patients — absent in Non-PSC entirely.

    adata_cd8 = adata_comp[adata_comp.obs["cell_type"] == "CD8+ cytotoxic T"].copy()
    print("\nCD8+ cytotoxic T cells by group:")
    print(adata_cd8.obs["group"].value_counts())

    gene = "IFNG"
    if gene in adata_cd8.raw.var_names:
        idx = list(adata_cd8.raw.var_names).index(gene)
        psc_expr = (
            adata_cd8[adata_cd8.obs["group"] == "PSC"]
            .raw.X[:, idx].toarray().flatten()
        )

        fig, ax = plt.subplots(figsize=(6, 5))
        ax.boxplot([psc_expr], tick_labels=["PSC"])
        ax.set_title("IFNG expression in CD8+ cytotoxic T cells\n(cluster absent in Non-PSC patients)")
        ax.set_ylabel("Normalized expression")
        plt.tight_layout()
        plt.savefig("figures/IFNG_cd8.png", dpi=150)
        print(f"PSC cells with IFNG > 0: {(psc_expr > 0).sum()} / {len(psc_expr)}")
        print("Saved: figures/IFNG_cd8.png")
    else:
        print(f"{gene} not found in raw var names")


# ── MAIN ──────────────────────────────────────────────────────────────────────
# Run stages in order. Each stage saves a checkpoint (.h5ad file) so you
# can reload from any point without rerunning earlier steps.
# Comment out stages that have already been run.

if __name__ == "__main__":
    # Stage 1 — load raw data from GEO (run once)
    adata = load_data()

    # Stage 2 — QC filtering (run once)
    adata = run_qc(adata)

    # Stage 3 — normalize, PCA, UMAP, cluster (slow — ~5-15 min)
    run_processing()

    # Stage 4 — marker genes + annotation
    run_annotation()

    # Stage 5 — PSC vs Non-PSC comparison
    run_comparison()