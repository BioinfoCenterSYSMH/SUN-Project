## Example script for SUNbeta12 using bundled HESC demo data
##
## This script:
## 1) loads bundled hesc_demo (counts + meta.data only),
## 2) applies UpdateSeuratObject for compatibility,
## 3) builds SNN graph with the project preprocessing sequence,
## 4) runs SUN on RNA_snn.

library(SUNbeta12)
library(Seurat)
data("hesc_demo", package = "SUNbeta12")
seu <- hesc_demo

# Some local Seurat setups may require object structure update.
seu <- tryCatch(
  SeuratObject::UpdateSeuratObject(seu),
  error = function(e) seu
)

set.seed(42)
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:10, verbose = FALSE)

if (!"RNA_snn" %in% names(seu@graphs)) {
  stop("RNA_snn graph not found. SUN requires a Seurat object with an SNN graph.", call. = FALSE)
}

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
