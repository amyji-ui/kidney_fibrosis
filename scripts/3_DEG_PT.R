# Qinglin Kong
# 3/13/26

# Load processed Seurat object
# Find cluster marker genes
# Subclustering analysis of PT in the UUO kidney

library(Seurat)
library(ggplot2)
library(harmony)

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("presto", quietly = TRUE))
  devtools::install_github("immunogenomics/presto")

library(presto)
library(dplyr)

out_obj_dir <- "results/seurat_objects"
out_fig_dir <- "results/figures"

ckd_pt <- readRDS("results/seurat_objects/02_gse182256_PT_subset.rds")
table(ckd_pt$orig.ident)

# figure 2 uses UUO kidneys
ckd_pt_uuo <- subset(ckd_pt, subset = orig.ident %in% c("UUO1", "UUO2"))
ckd_pt_uuo

# -----------------
# Run the pipeline again on UUO kieneys
# -----------------

ckd_pt_uuo <- NormalizeData(ckd_pt_uuo)

ckd_pt_uuo <- FindVariableFeatures(ckd_pt_uuo,
                                   selection.method = "vst",
                                   nfeatures = 2000)
top10 <- head(VariableFeatures(ckd_pt_uuo), 10)
p_var <- VariableFeaturePlot(ckd_pt_uuo)
p_var_lab <- LabelPoints(plot = p_var, points = top10, repel = TRUE)
p_var_lab
ggsave(
  filename = file.path(out_fig_dir, "03_PT_find_variable_features.png"),
  plot = p_var_lab,
  width = 6,
  height = 6,
  dpi = 300
)

ckd_pt_uuo <- ScaleData(ckd_pt_uuo,
                        vars.to.regress = c("nCount_RNA", "percent.mt"))

ckd_pt_uuo <- RunPCA(
  ckd_pt_uuo,
  assay = "RNA",
  npcs = 8,
  features = VariableFeatures(object = ckd_pt_uuo),
  ndims.print = 1:5,
  nfeatures.print = 10
)

# check where the elbow is to confirm if 10 dims is reasonable
p_elbow <- ElbowPlot(ckd_pt_uuo, ndims = 8)
p_elbow
ggsave(
  filename = file.path(out_fig_dir, "03_PT_pca_elbowplot.png"),
  plot = p_elbow,
  width = 6,
  height = 6,
  dpi = 300
)

# which genes contribute most to each PC
p_dim_loadings <- VizDimLoadings(ckd_pt_uuo, dims = 1:2, reduction = "pca")
p_dim_loadings
ggsave(
  filename = file.path(out_fig_dir, "03_PT_pca_dim_loadings.png"),
  plot = p_dim_loadings,
  width = 6,
  height = 6,
  dpi = 300
)

# show cells projected onto the first two PCs
p_dim_plot <- DimPlot(ckd_pt_uuo, reduction = "pca") + NoLegend()
p_dim_plot
ggsave(
  filename = file.path(out_fig_dir, "03_PT_pca_dim_plot.png"),
  plot = p_dim_plot,
  width = 6,
  height = 6,
  dpi = 300
)

ckd_pt_uuo <- RunHarmony(ckd_pt_uuo, group.by.vars = "orig.ident")

ckd_pt_uuo <- RunUMAP(ckd_pt_uuo, reduction = "harmony", dims = 1:8)
ckd_pt_uuo <- FindNeighbors(ckd_pt_uuo, reduction = "harmony", dims = 1:8)
ckd_pt_uuo <- FindClusters(ckd_pt_uuo, algorithm = "leiden", resolution = 0.5)

p_umap_clusters <- DimPlot(ckd_pt_uuo,
                           reduction = "umap",
                           label = TRUE,
                           pt.size = 0.5) + NoLegend()
p_umap_clusters
ggsave(
  filename = file.path(out_fig_dir, "03_PT_UUO_umap_seurat_clusters.png"),
  plot = p_umap_clusters,
  width = 6,
  height = 6,
  dpi = 300
)

# -----------------
# (under construction)
# Find differentially expressed features (cluster markers)
# -----------------

# find all differentially expressed markers
ckd_pt_markers <- FindAllMarkers(ckd_pt, 
                                 only.pos = TRUE, 
                                 logfc.threshold = 0.25)
ckd_pt_top10_markers <- ckd_pt.markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 10)


# save the markers to file
dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)

write.table(
  ckd_pt.markers,
  file = "results/markers/03_PT_all_cluster_markers.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ckd_pt.top10_markers,
  file = "results/markers/03_PT_top10_markers_per_cluster.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)