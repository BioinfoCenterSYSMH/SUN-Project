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

# seu: a Seurat object with a precomputed SNN graph (e.g., RNA_snn)
res <- SUN(
  seurat_obj = seu,
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

- Input must be a Seurat object with an available SNN graph.
- `SUN_Label` is written to `meta.data` for downstream analysis and visualization.

## License

MIT
