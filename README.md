# CAR-T Cytopenia: scRNA-seq Analysis of Bone Marrow

Independent reanalysis of bone marrow single-cell RNA sequencing data from patients who developed prolonged cytopenia following CD19 CAR-T cell therapy for relapsed/refractory large B-cell lymphoma (rrLBCL).

## Background

CAR-T cell therapy can eliminate cancer that has failed all other treatments. But in approximately 30% of patients, the bone marrow fails to recover afterward — a complication called prolonged cytopenia (PC). White blood cell counts stay dangerously low for months. Patients die from infections, not from cancer.

The biology of why some patients' bone marrow recovers and others' doesn't is an active and unsolved research problem.

## Question

What is different at the single-cell level in the bone marrow of patients whose marrow fails to recover after CAR-T therapy?

## Dataset

**GEO Accession: GSE216005**

Single-cell RNA sequencing of bone marrow aspirates from 22 individuals:
- 11 rrLBCL patients with prolonged cytopenia (PSC) after CAR-T
- 5 rrLBCL patients without prolonged cytopenia (Non-PSC)
- 5 healthy bone marrow donors
- 1 chemo-associated cytopenia control
- 2 excluded (bone marrow involvement with AML/LBCL)

Total: 106,198 cells after QC filtering across 27,089 genes.

**Source paper:** Strati P, Li X, Deng Q et al. *Prolonged cytopenia following CD19 CAR T cell therapy is linked with bone marrow infiltration of clonally expanded IFNγ-expressing CD8 T cells.* Cell Reports Medicine. 2023;4(8):101158. PMID: 37586321

## Key Finding

CD8+ cytotoxic T cells are massively expanded in PSC patients compared to Non-PSC patients (14.9% vs 0% of bone marrow cells). This population — characterized by high expression of NKG7, CD8A, GZMH, and CCL5 — is entirely absent in patients whose bone marrow recovered.

Within this PSC-exclusive population, 11.5% of cells have detectable IFNG (IFN-γ) expression. IFN-γ suppresses hematopoietic stem cell self-renewal and differentiation, providing a mechanistic explanation for prolonged bone marrow failure.

This reproduces the core finding of Strati et al. 2023 using an independent pipeline.

## Results

| Figure | Description |
|--------|-------------|
| `figures/umap_clusters.png` | 25 Leiden clusters from unsupervised clustering |
| `figures/umap_annotated.png` | Cell type annotations across all clusters |
| `figures/umap_groups.png` | Patient group distribution across UMAP |
| `figures/umap_groups_annotated.png` | PSC vs Non-PSC vs Healthy overlay |
| `figures/proportion_comparison.png` | Cell type enrichment: PSC vs Non-PSC |
| `figures/IFNG_cd8.png` | IFN-γ expression in CD8+ cytotoxic T cells |

## Pipeline

```
Raw 10x count matrices (GEO: GSE216005)
    ↓
QC filtering (min 200 genes, max 15% mito, max 4000 genes)
    ↓
Normalization (10k counts/cell, log1p)
    ↓
Highly variable genes (batch_key="patient")
    ↓
PCA (50 components) → Neighbors → UMAP
    ↓
Leiden clustering (resolution 0.4 → 25 clusters)
    ↓
Marker gene annotation
    ↓
PSC vs Non-PSC proportion analysis + IFNG expression
```

## Repository Structure

```
├── analysis/
│   └── cart_analysis.py     # main analysis — 5 stages, checkpoint saves
├── notebooks/
│   └── pbmc3k.py            # tutorial pipeline (PBMC 3k dataset)
├── figures/                 # all output plots
└── README.md
```

## How to Run

**Requirements:**
```bash
conda create -n scenv python=3.11
conda activate scenv
pip install scanpy
conda install -c conda-forge python-igraph
```

**Data:**
Download GSE216005 from NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216005

Unpack into `data/GSE216005/` so each patient folder (`pt01/` through `pt22/`) contains `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`.

Update `DATA_DIR` in `cart_analysis.py` to match your local path.

**Run:**
```bash
python analysis/cart_analysis.py
```

Each stage saves a checkpoint `.h5ad` file. To resume from a checkpoint, comment out earlier stages in the `__main__` block.

## Context

This analysis was motivated by a personal interest in understanding why CAR-T therapy — a treatment that can eliminate cancer — sometimes fails patients at the final step. The bone marrow failure that kills some patients after successful CAR-T treatment is preventable. Understanding it at the molecular level is the first step toward preventing it.

## Author

Joshua Brilon  
B.S. Biology (anticipated June 2028), University of Oregon  
Minors: Bioengineering, Computer Science  
GitHub: [github.com/JoshuaBrilon](https://github.com/JoshuaBrilon)
