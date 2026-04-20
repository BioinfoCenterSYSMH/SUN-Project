## SUNbeta12 v0.1.2

### Highlights
- Added bundled **HESC demo dataset** (`hesc_demo`) for out-of-the-box examples.
- Updated the example to rebuild the SNN graph from counts + metadata before running SUN.
- Clarified that `SUN()` requires a Seurat object with an available SNN graph (e.g., `RNA_snn`).
- Added compatibility handling with `SeuratObject::UpdateSeuratObject()` in the demo workflow.

### HESC Demo Workflow
The bundled HESC example follows the project preprocessing sequence:
1. `NormalizeData`
2. `FindVariableFeatures`
3. `ScaleData`
4. `RunPCA`
5. `FindNeighbors`
6. `SUN(...)`

### Packaging
- Bumped package version to **0.1.2**
- Source archive: `SUNbeta12_0.1.2.tar.gz`

### Installation
```r
install.packages("remotes")
remotes::install_github("BioinfoCenterSYSMH/SUN-Project@v0.1.2")
```
