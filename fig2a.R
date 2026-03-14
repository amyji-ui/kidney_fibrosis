library(Matrix)
library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)

# create the Seurat object from the counts export 
mat_path  <- "data/GSE182256_Export_counts.rds"
mat <- readRDS(mat_path)
ckd <- CreateSeuratObject(counts = mat, project = "kidney_fibrosis", min.cells = 3, min.features = 200)  # nolint: line_length_linter.
print(ckd)

# preprocessing
ckd[["nFeature_RNA"]] <- colSums(GetAssayData(ckd, layer = "counts") > 0)
ckd[["nCount_RNA"]] <- colSums(GetAssayData(ckd, layer = "counts"))
ckd[["percent.mt"]] <- PercentageFeatureSet(ckd, pattern = "^MT-")

# visualize QC metrics
vp <- VlnPlot(ckd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) # nolint: line_length_linter.
ggsave("output/qc_violin.png", plot = vp, width = 8, height = 6, dpi = 300)