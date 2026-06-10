# Fetal Uterus Development Atlas

Code and analysis pipeline for the fetal uterus development single-cell atlas. This repository contains the core analysis scripts used for snRNA-seq cell-type clustering, MNN batch integration, spatial transcriptomics (STAGATE), WGCNA co-expression network analysis, pseudotime trajectory inference, and cell-cell communication (CellChat).

---

## 1. System Requirements

### Operating System

| Item | Version / Detail |
|---|---|
| OS | Linux (CentOS 7, kernel 3.10.0-1160.el7.x86\_64) |
| R | 4.2.0 (2022-04-22) |
| Python | 3.12.4 |

### R Package Dependencies

| Package | Version | Usage |
|---|---|---|
| Seurat | 5.0.1 | Core single-cell analysis (normalization, PCA, clustering, UMAP, MNN integration) |
| SeuratObject | 5.0.1 | Seurat data structure (auto-loaded with Seurat) |
| SeuratWrappers | 0.3.2 | Interface wrappers for Monocle3, FastMNN etc. |
| Matrix | 1.6.4 | Sparse matrix operations |
| dplyr | 1.1.4 | Data manipulation |
| ggplot2 | 3.5.1 | Visualization |
| patchwork | 1.2.0 | Multi-panel figure composition |
| RColorBrewer | 1.1.3 | Color palettes |
| hdWGCNA | 0.4.0 | High-dimensional WGCNA co-expression network analysis |
| UCell | 2.10.1 | Module scoring for WGCNA |
| SCENIC | 1.3.0 | Gene regulatory network inference |
| CellChat | 2.1.2 | Cell-cell communication analysis |
| CytoTRACE2 | 1.1.0 | Differentiation potential estimation |
| monocle3 | 1.4.26 | Pseudotime trajectory inference |
| rbokeh | 0.5.2 | Interactive visualization |
| nlme | 3.1.163 | Non-linear mixed effects models |
| curl | 5.2.3 | Network data transfer |
| cluster | 2.1.6 | Clustering algorithms |
| stringi | 1.8.3 | String processing |
| stringr | 1.5.1 | String manipulation |
| future | 1.33.1 | Parallel processing backend |

### Python Package Dependencies

(Required for Coding-3: STAGATE spatial analysis, and Coding-4: spatial visualization)

| Package | Version | Usage |
|---|---|---|
| Python | 3.7.12 / 3.12.4 | |
| scanpy | 1.9.3 / 1.11.5 | Single-cell analysis in Python |
| anndata | 0.8.0 / 0.10.9 | Annotated data structure |
| pandas | 2.2.2 | Data manipulation |
| numpy | 1.26.4 | Numerical computation |
| matplotlib | 3.8.4 | Visualization |
| scipy | 1.13.1 | Scientific computing |
| Pillow | 10.3.0 | Image I/O |
| STAGATE | 1.0.1 | Spatial transcriptomics graph attention network |
| tensorflow | 1.15.0 | Deep learning backend for STAGATE |
| opencv-python | 4.8.0 | Cell boundary detection (Coding-4) |


### Non-standard Hardware

- A workstation/server with ≥64 GB RAM is recommended for processing the full snRNA-seq dataset (61,099 genes × 18,718 cells for the EPI subset; larger for the merged fetal+adult dataset).
- A CUDA-capable GPU may accelerate STAGATE training (optional).

---

## 2. Installation Guide

### R Packages

```r
# Core packages from CRAN
install.packages(c("Seurat", "Matrix", "curl", "nlme", "dplyr", "ggplot2", 
                   "RColorBrewer", "patchwork", "rbokeh", "stringi", "stringr", "future"))

# Bioconductor / GitHub packages
# SeuratWrappers, hdWGCNA, SCENIC, CellChat, UCell, monocle3, CytoTRACE2
```

Typical install time: ~30–60 minutes on a standard desktop.

### Python Packages

STAGATE (Coding-3) requires Python 3.7 and tensorflow 1.15.0:

```bash
pip install scanpy anndata pandas numpy matplotlib scipy pillow
pip install stagate-pyg tensorflow==1.15.0
```

Spatial visualization (Coding-4) additionally requires:

```bash
pip install opencv-python
```

Typical install time: ~20–40 minutes.

---

## 3. Demo

### Demo Data

A demo Seurat object is provided at [`EPI.rds`](https://github.com/dingli-pku/Fetal-uterus-development-atlas/blob/main/EPI.rds) (547 MB, epithelial cell subset).  
It contains 61,099 genes × 18,718 cells with the following pre-loaded components:

| Slot | Content |
|---|---|
| `RNA` assay | Log-normalized expression data |
| `mnn.reconstructed` assay | MNN-reconstructed expression |
| `pca` reduction | 50 PCs |
| `integrated.mnn` reduction | FastMNN-corrected embeddings (50 dims) |
| `umap.mnn` reduction | 2D UMAP projection |
| `cell_sub_type` metadata | Cell subtype annotation |

### Demo Execution

The demo applies the **general clustering pipeline** (Coding-1) to the EPI subset with **EPI-specific** parameters.

```r
# ---- Environment ----
library(Seurat)
library(SeuratWrappers)

# ---- Load demo data ----
EPI <- readRDS("EPI.rds")

# ---- Parameter reference (from EPI@commands) ----
# These are the actual parameters used, recorded in the Seurat object:

# Step 1: NormalizeData
EPI <- NormalizeData(EPI, assay = "RNA")
#   method: LogNormalize, scale.factor: 10000

# Step 2: FindVariableFeatures
EPI <- FindVariableFeatures(EPI, assay = "RNA", selection.method = "vst", nfeatures = 2000)

# Step 3: ScaleData
EPI <- ScaleData(EPI, assay = "RNA")
#   model.use: linear, scale.max: 10

# Step 4: Run PCA
EPI <- RunPCA(EPI, assay = "RNA", npcs = 50, seed.use = 42)

# Step 5: FindNeighbors on PCA space — dims = 1:10 (EPI-specific)
EPI <- FindNeighbors(EPI, reduction = "pca", dims = 1:10)

# Step 6: FastMNN integration
EPI <- IntegrateLayers(
  object = EPI,
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

# Step 7: FindNeighbors on MNN space — dims = 1:8 (EPI-specific)
EPI <- FindNeighbors(EPI, reduction = "integrated.mnn", dims = 1:8)

# Step 8: Run UMAP on MNN space — dims = 1:8, metric = cosine (EPI-specific)
EPI <- RunUMAP(
  EPI,
  reduction      = "integrated.mnn",
  dims           = 1:8,
  n.neighbors    = 30,
  min.dist       = 0.3,
  metric         = "cosine",
  seed.use       = 42,
  reduction.name = "umap.mnn"
)

# Step 9: Clustering — resolution = 2 (EPI-specific)
EPI <- FindClusters(EPI, resolution = 2, cluster.name = "mnn_clusters")

# Step 10: Visualization
DimPlot(EPI, reduction = "umap.mnn", group.by = "cell_sub_type", label = TRUE)
```

### Expected Output

- A UMAP plot colored by `cell_sub_type`.

### Expected Runtime

Approximately **3–5 minutes** on a standard desktop/laptop with 16 GB RAM for the full pipeline on the EPI subset (18,718 cells).

---

## 4. Instructions for Use

### How to Run on Your Own Data

**Prerequisite**: Your input should be a **Seurat object** (`.rds` or `.RData`) containing raw or normalized RNA counts.

1. **Place your Seurat object** in your working directory.
2. **Adjust parameters** (`dims`, `resolution`) as needed for your cell lineage.
3. **Source the pipeline code** (see Coding-1-snRNA_fetal.txt for the general template, or use the demo code above with your parameters).

### Analysis Overview

Each script in this repository addresses a specific analysis task:

| Script | Task | Environment |
|---|---|---|
| `Coding-1-snRNA_fetal.txt` | snRNA-seq fetal cell-type clustering + MNN integration | R (Seurat 5) |
| `Coding-2-snRNA_fetal_adult.txt` | Joint fetal + adult MNN integration | R (Seurat 5) |
| `Coding-3-BMK-STAGATE.txt` | Spatial transcriptomics (STAGATE) | Python (scanpy / STAGATE) |
| `Coding-4-BMK-visualization.txt` | Spatial cell-type boundary visualization | Python (opencv / scipy) |
| `Coding-5-WGCNA.txt` | Co-expression network analysis (hdWGCNA) | R (hdWGCNA) |
| `Coding-6-clustering.txt` | Correlation-based hierarchical clustering | R |
| `Coding-7-Pseudotime.txt` | Pseudotime trajectory inference (monocle3) | R (monocle3) |
| `Coding-8-Pseudotime_gene.txt` | Pseudotime gene dynamics (GAM smoothing) | R (ggplot2) |
| `Coding-9-Potential.txt` | Differentiation potential (CytoTRACE2) | R (CytoTRACE2) |
| `Coding-10-CellChat.txt` | Cell-cell communication (CellChat) | R (CellChat) |

### Reproduction of Results

To reproduce results from the manuscript:
1. Run Coding-1 or Coding-2 to obtain clustered, MNN-integrated Seurat objects.
2. Subset each major cell type (EPI, stromal, etc.) and re-run the clustering pipeline with lineage-specific parameters.
3. Proceed with downstream analyses: WGCNA (Coding-5), pseudotime (Coding-7/8/9), CellChat (Coding-10), and clustering (Coding-6).
4. For spatial data, run Coding-3 (STAGATE) then Coding-4 (visualization).

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details (Copyright © 2025 dingli-pku).

