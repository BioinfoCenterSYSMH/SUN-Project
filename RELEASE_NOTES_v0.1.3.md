# SUNbeta13 v0.1.3 (sun-beta-1.3)

## Highlights

- **Single implementation file** `R/SUN_beta13.R`: self-contained SUN beta 1.3 (no `source()` of other in-repo R scripts).
- **Recommendation scoring** uses **`bluster::approxSilhouette`** on a linear embedding from `Seurat::Embeddings()` (default reduction from env `SUN_SILHOUETTE_REDUCTION` or `"pca"`), not graph-shortest-path silhouette.
- **Plateau–cache-aligned** stable intervals and **`SUN()`** workflow aligned with the prior `SUN_beta1.1.R` + `SUN_beta_1.3.R` chain.

## Installation

Install **bluster** from Bioconductor, then:

```r
remotes::install_github("BioinfoCenterSYSMH/SUN-Project", ref = "sun-beta-1.3")
```

## Documentation

- See **README.md** on branch `sun-beta-1.3` for usage, including `silhouette_reduction` / `silhouette_dims` matching `FindNeighbors`.
