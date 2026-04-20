# SUNbeta12

`SUNbeta12` is an R package wrapper for the SUN beta 1.2 resolution-recommendation workflow on Seurat SNN graphs.

## What it does

- scans resolution intervals from `r_min` to `r_max`
- reuses and validates existing interval cache when available
- applies beta 1.2 plateau-cache aligned interval replacement
- recommends a resolution using graph-distance silhouette + interval width
- writes `SUN_Label` (and compatibility label `SUN_Label_beta11`) into `meta.data`

## Install (after pushing to GitHub)

```r
install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/SUNbeta12")
```

## Minimal usage

```r
library(SUNbeta12)
library(Seurat)

# seu should already contain an SNN graph, e.g. RNA_snn
res <- SUN(
  seurat_obj = seu,
  r_min = 0,
  r_max = 1,
  tools_key = "all_interval_analysis",
  graph.name = "RNA_snn",
  random.seed = 42L
)

seu_out <- res$seurat_obj
rec_r <- res$recommended_resolution
head(seu_out@meta.data$SUN_Label)
```

## Build source tarball

From terminal:

```bash
R CMD build SUNbeta12
```

This produces a `.tar.gz` package you can upload/release.
