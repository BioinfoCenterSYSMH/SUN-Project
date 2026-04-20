# SUN-Project

SUN provides reproducible resolution recommendation for Seurat SNN clustering using stable intervals and graph-distance silhouette.

## Overview

SUN is an R toolkit for automatic and auditable resolution selection in graph-based single-cell clustering workflows built on Seurat objects.

## Installation

```r
install.packages("remotes")
remotes::install_github("BioinfoCenterSYSMH/SUN-Project")
```

## Quick Start

```r
library(SUNbeta12)
library(Seurat)
data("hesc_demo", package = "SUNbeta12")
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

# Run SUN
res <- SUN(
  seurat_obj = seurat_obj,
  r_min = 0,
  r_max = 1,
  tools_key = "all_interval_analysis",
  graph.name = "RNA_snn",
  random.seed = 42L
)

# Recommended resolution
res$recommended_resolution

# Output Seurat object with SUN labels
seu_out <- res$seurat_obj
head(seu_out@meta.data$SUN_Label)
```

## Output

`SUN()` returns a list containing:
- `seurat_obj`: Seurat object with SUN results written back
- `recommended_resolution`: recommended Louvain resolution
- `recommendation`: ranking details for candidate intervals

## Notes

- `SUN()` accepts Seurat objects with an available SNN graph (e.g., `RNA_snn`).
- If no SNN graph exists, run preprocessing first. For the bundled HESC demo,
  the SNN construction steps are:
  `NormalizeData` -> `FindVariableFeatures` -> `ScaleData` -> `RunPCA` -> `FindNeighbors`.
- If object-version compatibility errors appear, apply
  `SeuratObject::UpdateSeuratObject()` before preprocessing.
- `SUN_Label` is written to `meta.data` for downstream analysis and visualization.

## License

MIT
