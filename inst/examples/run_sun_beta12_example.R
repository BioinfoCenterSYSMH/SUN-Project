## Example script for SUNbeta12
##
## Prerequisite:
## - `seu` is a Seurat object that already has an SNN graph
##   (e.g., `RNA_snn` after FindNeighbors)

library(SUNbeta12)

out <- SUN(
  seurat_obj = seu,
  r_min = 0,
  r_max = 1,
  tools_key = "all_interval_analysis",
  graph.name = "RNA_snn",
  random.seed = 42L
)

seu_out <- out$seurat_obj
print(out$recommended_resolution)
table(seu_out@meta.data$SUN_Label, useNA = "ifany")
