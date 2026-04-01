# Qinglin Kong
# 3/12/26

# Set paths
# Load counts and create Seurat object
# Load and attach reference metadata
# Compute QC metrics
# Make QC plots and inspect deposited filtering
# Use deposited QC-filtered dataset for downstream analysis

# Normalize data
# Find highly variable features
# Scale data
# Run PCA

# Run Harmony integration

# Run UMAP and clustering
# Save processed Seurat object

library(Seurat)
library(ggplot2)
library(harmony)

# -----------------
# Set paths
# -----------------

mat_path  <- "data/GSE182256_Export_counts.rds"
meta_path <- "data/GSE182256_Export_Metadata.txt"
umap_path <- "data/GSE182256_EXPORT_umap.txt"

out_obj_dir <- "results/seurat_objects"
out_fig_dir <- "results/figures"

dir.create(out_obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------
# Load counts and create object
# -----------------

mat <- readRDS(mat_path)
ckd <- CreateSeuratObject(counts = mat, 
                          project = "GSE182256", 
                          min.cells = 3, 
                          min.features = 200)
print(ckd)

# -----------------
# Load and attach metadata to Seurat obj
# -----------------
meta <- read.delim(meta_path, 
                   header = TRUE, 
                   row.names = 1, 
                   sep = "\t", 
                   check.names = FALSE)
meta <- meta[, c("orig.ident", "Cluster"), drop = FALSE]
colnames(meta)[colnames(meta) == "Cluster"] <- "paper_cluster"

common_cells <- intersect(colnames(ckd), rownames(meta))
ckd <- subset(ckd, cells = common_cells)
meta <- meta[common_cells, , drop = FALSE]
ckd <- AddMetaData(ckd, metadata = meta)

head(ckd@meta.data)
table(ckd$paper_cluster)

# -----------------
# Compute QC metrics
# -----------------
ckd[["percent.mt"]]   <- PercentageFeatureSet(ckd, pattern = "^mt-")
head(ckd@meta.data)

# -----------------
# Make QC plots
# -----------------
# violin plots
p_vln <- VlnPlot(ckd, 
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 ncol = 3, 
                 pt.size = 0)
p_vln
ggsave(
  filename = file.path(out_fig_dir, "01_qc_violin_before_filtering.png"),
  plot = p_vln,
  width = 9,
  height = 6,
  dpi = 300
)

# check relationship of UMI counts vs number of genes
p_scatter1 <- FeatureScatter(ckd, 
                             feature1 = "nCount_RNA", 
                             feature2 = "nFeature_RNA")
p_scatter1
ggsave(
  filename = file.path(out_fig_dir, "01_qc_scatter_counts_vs_features.png"),
  plot = p_scatter1,
  width = 6,
  height = 6,
  dpi = 300
)

# check if high mitochondrial cells also have low counts
p_scatter2 <- FeatureScatter(ckd, 
                             feature1 = "nCount_RNA", 
                             feature2 = "percent.mt")
p_scatter2
ggsave(
  filename = file.path(out_fig_dir, "01_qc_scatter_counts_vs_percent_mt.png"),
  plot = p_scatter2,
  width = 6,
  height = 6,
  dpi = 300
)

# summary(ckd$nFeature_RNA)
# summary(ckd$percent.mt)

# -----------------
# Use deposited QC-filtered dataset for downstream analysis
# -----------------

# ckd_filt <- subset(ckd, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 50)
ckd_filt <- ckd

# -----------------
# Normalize data
# -----------------

ckd_filt <- NormalizeData(ckd_filt)

# -----------------
# Find highly variable features
# -----------------

ckd_filt <- FindVariableFeatures(ckd_filt, 
                                 selection.method = "vst", 
                                 nfeatures = 2000)

top10 <- head(VariableFeatures(ckd_filt), 10)
p_var <- VariableFeaturePlot(ckd_filt)
p_var_lab <- LabelPoints(plot = p_var, points = top10, repel = TRUE)
p_var_lab
ggsave(
  filename = file.path(out_fig_dir, "01_find_variable_features.png"),
  plot = p_var_lab,
  width = 6,
  height = 6,
  dpi = 300
)

# -----------------
# Scale data
# -----------------

ckd_filt <- ScaleData(ckd_filt, vars.to.regress = c("nCount_RNA", "percent.mt"))

# -----------------
# Run PCA
# -----------------

# using the paper's code
ckd_filt <- RunPCA(
  ckd_filt, 
  assay = "RNA", 
  npcs = 30, 
  features = VariableFeatures(object = ckd_filt), 
  ndims.print = 1:5, 
  nfeatures.print = 10
)

# check where the elbow is to confirm if 30 dims is reasonable
p_elbow <- ElbowPlot(ckd_filt, ndims = 30)
ggsave(
  filename = file.path(out_fig_dir, "01_pca_elbowplot.png"),
  plot = p_elbow,
  width = 6,
  height = 6,
  dpi = 300
)
p_elbow

# which genes contribute most to each PC
p_dim_loadings <- VizDimLoadings(ckd_filt, dims = 1:2, reduction = "pca")
p_dim_loadings
ggsave(
  filename = file.path(out_fig_dir, "01_pca_dim_loadings.png"),
  plot = p_dim_loadings,
  width = 6,
  height = 6,
  dpi = 300
)

# show cells projected onto the first two PCs
p_dim_plot <- DimPlot(ckd_filt, reduction = "pca") + NoLegend()
p_dim_plot
ggsave(
  filename = file.path(out_fig_dir, "01_pca_dim_plot.png"),
  plot = p_dim_plot,
  width = 6,
  height = 6,
  dpi = 300
)

# -----------------
# Run Harmony integration
# -----------------

ckd_filt <- RunHarmony(ckd_filt, group.by.vars = "orig.ident")

# -----------------
# Run UMAP and clustering
# -----------------

ckd_filt <- RunUMAP(ckd_filt, reduction = "harmony", dims = 1:30)
ckd_filt <- FindNeighbors(ckd_filt, reduction = "harmony", dims = 1:30)
ckd_filt <- FindClusters(ckd_filt)

p_umap_clusters <- DimPlot(ckd_filt, 
                           reduction = "umap", 
                           label = TRUE, 
                           pt.size = 0.5) + NoLegend()
p_umap_clusters
ggsave(
  filename = file.path(out_fig_dir, "01_umap_seurat_clusters.png"),
  plot = p_umap_clusters,
  width = 6,
  height = 6,
  dpi = 300
)

# -----------------
# Save processed Seurat object
# -----------------
saveRDS(
  ckd_filt,
  file = file.path(out_obj_dir, "01_umap_gse182256_processed.rds")
)

# save cluster numbers and UMAP coordinates as a table
umap_df <- as.data.frame(Embeddings(ckd_filt, "umap"))
umap_df$cell <- rownames(umap_df)
umap_df$seurat_clusters <- ckd_filt$seurat_clusters
umap_df$paper_cluster <- ckd_filt$paper_cluster
umap_df$orig.ident <- ckd_filt$orig.ident

write.table(
  umap_df,
  file = file.path(out_obj_dir, "01_umap_coordinates_and_clusters.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
