# SUN-Project

SUN provides reproducible resolution recommendation for Seurat SNN clustering using stable intervals and graph-distance silhouette.

## Overview
SUN is an R toolkit for automatic and auditable resolution selection in graph-based single-cell clustering workflows.

## Key Features
- Stable-interval discovery across resolution scans
- Graph-distance silhouette scoring on SNN-induced subgraphs
- Deterministic recommendation and cache-aware workflow
- Seamless integration with Seurat objects

## Installation
```r
install.packages("remotes")
remotes::install_github("BioinfoCenterSYSMH/SUN-Project")
library(SUNbeta12)

res <- SUN(
  seurat_obj = seu,
  r_min = 0,
  r_max = 1,
  tools_key = "all_interval_analysis",
  graph.name = "RNA_snn",
  random.seed = 42L
)

res$recommended_resolution
head(res$seurat_obj@meta.data$SUN_Label)
