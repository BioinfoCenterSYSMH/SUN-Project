# SUN-Project (SUN beta 1.3)

SUN provides reproducible resolution recommendation for Seurat SNN clustering using **stable resolution intervals**, **plateau–cache-aligned** interval replacement, and **embedding-based silhouette** scoring via **`bluster::approxSilhouette`** (beta 1.3). This replaces the older graph-distance silhouette used in earlier beta lines for the ranking step.

Branch: [`sun-beta-1.3`](https://github.com/BioinfoCenterSYSMH/SUN-Project/tree/sun-beta-1.3).

## Overview

SUN is an R toolkit for automatic and auditable resolution selection in graph-based single-cell clustering workflows built on Seurat objects. The implementation on this branch lives primarily in **`R/SUN_beta13.R`** (standalone: no `source()` of other project R files).

## Installation

Bioconductor **bluster** is required:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bluster")

install.packages("remotes")
remotes::install_github("BioinfoCenterSYSMH/SUN-Project", ref = "sun-beta-1.3")
```

## Quick start

```r
library(SUNbeta13)
library(Seurat)
data("hesc_demo", package = "SUNbeta13")
seurat_obj <- hesc_demo

# Compatibility guard for local Seurat object structure
seurat_obj <- tryCatch(
  SeuratObject::UpdateSeuratObject(seurat_obj),
  error = function(e) seurat_obj
)

# Build SNN graph (same logic as the project workflow)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:10, verbose = FALSE)

# Run SUN (beta 1.3): match silhouette embedding/dims to FindNeighbors
res <- SUN(
  seurat_obj = seurat_obj,
  r_min = 0,
  r_max = 1,
  tools_key = "all_interval_analysis",
  graph.name = "RNA_snn",
  random.seed = 42L,
  silhouette_reduction = "pca",
  silhouette_dims = 1:10
)

res$recommended_resolution

seu_out <- res$seurat_obj
head(seu_out@meta.data$SUN_Label)
```

## Output

`SUN()` returns a list containing:

- `seurat_obj`: Seurat object with SUN results written back
- `recommended_resolution`: recommended Louvain resolution
- `recommendation`: ranking details for candidate intervals (embedding silhouette in beta 1.3)

`obj@tools$SUN` stores summary metadata (including `silhouette_reduction` / `silhouette_dims` when set).

## Notes

- `SUN()` expects an SNN graph (e.g. `RNA_snn`) and a reduction used for silhouette (default follows env `SUN_SILHOUETTE_REDUCTION` or `"pca"`). Use **`silhouette_reduction` / `silhouette_dims` consistent with `FindNeighbors(..., reduction=, dims=)`** (e.g. Harmony / MNN integrated embeddings).
- Legacy scripts `R/SUN_beta11.R` and `R/SUN_beta12.R` are kept in the repository for reference; **the installable package code path should rely on `R/SUN_beta13.R`** (see `.Rbuildignore` if present).
- If object-version compatibility errors appear, apply `SeuratObject::UpdateSeuratObject()` before preprocessing.
- `SUN_Label` and `SUN_Label_beta11` are written to `meta.data` for downstream analysis.

## License

MIT
