# Qinglin Kong
# 3/13/26

# Load processed Seurat object
# Find cluster marker genes
# Subclustering analysis of PT in the UUO kidney
# Find differentially expressed genes (cluster markers) in each cluster
# Find which clusters represent profibrotic PT cells.

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
out_mark_dir <- "results/markers"


# -----------------
# Helper functions
# -----------------

save_plot <- function(plot_obj, filename, width = 6, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(out_fig_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

save_table <- function(df, filename) {
  write.table(
    df,
    file = file.path(out_mark_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

make_bubble_plot <- function(obj, features, group.by = "seurat_clusters",
                             low_color = "white", high_color = "deeppink4",
                             dot.min = 0.05) {
  features <- features[features %in% rownames(obj)]
  DotPlot(
    obj,
    features = features,
    group.by = group.by,
    cols = c(low_color, high_color),
    dot.min = dot.min
  ) + RotatedAxis()
}

find_top_markers_per_cluster <- function(marker_table, top_n = 5) {
  marker_table |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = top_n)
}

clean_marker_table_for_plotting <- function(marker_table) {
  marker_table |>
    filter(!grepl("^mt-", gene, ignore.case = TRUE)) |>
    filter(!grepl("^Gm", gene)) |>
    filter(!grepl("Rik$", gene)) |>
    filter(!gene %in% c("Xist", "Tsix"))
}

# -----------------
# Load processed Seurat object
# -----------------

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
  filename = file.path(out_fig_dir, "03_PT_UUO_find_variable_features.png"),
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
  npcs = 30,
  features = VariableFeatures(object = ckd_pt_uuo),
  ndims.print = 1:5,
  nfeatures.print = 10
)

# check where the elbow is
p_elbow <- ElbowPlot(ckd_pt_uuo, ndims = 30)
p_elbow
ggsave(
  filename = file.path(out_fig_dir, "03_PT_UUO_pca_elbowplot.png"),
  plot = p_elbow,
  width = 6,
  height = 6,
  dpi = 300
)

# which genes contribute most to each PC
p_dim_loadings <- VizDimLoadings(ckd_pt_uuo, dims = 1:2, reduction = "pca")
p_dim_loadings
ggsave(
  filename = file.path(out_fig_dir, "03_PT_UUO_pca_dim_loadings.png"),
  plot = p_dim_loadings,
  width = 6,
  height = 6,
  dpi = 300
)

# show cells projected onto the first two PCs
p_dim_plot <- DimPlot(ckd_pt_uuo, reduction = "pca") + NoLegend()
p_dim_plot
ggsave(
  filename = file.path(out_fig_dir, "03_PT_UUO_pca_dim_plot.png"),
  plot = p_dim_plot,
  width = 6,
  height = 6,
  dpi = 300
)

ckd_pt_uuo <- RunHarmony(ckd_pt_uuo, group.by.vars = "orig.ident")

ckd_pt_uuo <- RunUMAP(ckd_pt_uuo, reduction = "harmony", dims = 1:30)
ckd_pt_uuo <- FindNeighbors(ckd_pt_uuo, reduction = "harmony", dims = 1:30)
ckd_pt_uuo <- FindClusters(ckd_pt_uuo)

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
# Bubble plot of representative genes in each cluster.
# -----------------

pt_marker_genes_fig_2b <- c(
  "Slc34a1",  # PT
  "Slc13a3",  # PT
  "Apom",     # precursor
  "Slc5a12",  # S1
  "Gsta4",    # S2
  "Slc6a13",  # S3
  "Pdgfb",    # profibrotic PT
  "Wfdc15b",  # transient mix
  "Mki67",    # proliferating
  "Cd52"      # immune
)

p_pt_marker_genes_bubble <- DotPlot(
  ckd_pt_uuo,
  features = pt_marker_genes_fig_2b,
  group.by = "seurat_clusters",
  cols = c("white", "deeppink4"),
  dot.min = 0.1
) + RotatedAxis()
p_pt_marker_genes_bubble
ggsave(
  filename = file.path(out_fig_dir, "03_PT_UUO_fig_2b_bubble_plot.png"),
  plot = p_pt_marker_genes_bubble,
  width = 6,
  height = 6,
  dpi = 300
)

# -----------------
# Find which clusters may represent profibrotic PT cells.
# -----------------
pt_profibrotic_panel <- c(
  "Lrp2", "Slc34a1", # PT identity
  "Pdgfb",           # profibrotic driver
  "Cd74", "Tnfrsf12a", "Cxcl1", "Cxcl10", "Cxcl16" # proinflammatory sig
)
p_profibrotic_bubble <- make_bubble_plot(
  ckd_pt_uuo,
  features = pt_profibrotic_panel
)
p_profibrotic_bubble
save_plot(p_profibrotic_bubble, "03_PT_UUO_profibrotic_panel.png")
# conclusion: cluster 10 might be profibrotic PT.

# -----------------
# Find which clusters may represent profibrotic tubule.
# -----------------
profib_tubule_genes <- c(
  "Pdgfb",    # profibrotic
  "Cxcl1",    # inflammatory
  "Epcam",    # epithelial / tubule
  "Lrp2",     # weak PT identity
  "Hnf4a",    # weak PT identity
  "Slc34a1",  # weak PT identity
  "Slc13a3",  # expected to be low / absent
  "Slc22a6",  # expected to be low / absent
  "Aqp1",     # descending loop / tubule-adjacent
  "Cryab",    # descending loop / stress
  "Cp"        # descending loop marker mentioned in paper
)
# keep only genes that actually exist in the Seurat object
profib_tubule_genes <- profib_tubule_genes[profib_tubule_genes %in% 
                                             rownames(ckd_pt_uuo)]

p_profib_tubule_bub <- make_bubble_plot(
  ckd_pt_uuo,
  features = profib_tubule_genes,
  dot.min = 0.1
)
p_profib_tubule_bub
save_plot(p_profib_tubule_bub, "03_PT_UUO_profibrotic_tubule_bubble_plot.png")

# conclusion: don't really see a clear cluster that represents proximal tubules 
# as described in the paper.

# ----------------
# Identify the cluster that are profibrotic PTs. 
# ----------------

# check gene enrichments, genes from two papers. 
# 1. Dove et al. 2022
# 2. Kirita et al. 2020
pt_marker_sets <- list(
  Precursor      = c("Apom"),
  S1             = c("Slc5a12", "Adra1a","Slc2a2", "Cubn"),
  S2             = c("Gsta4", "Slc22a30", "Slc7a13"),
  S3             = c("Slc6a13", "Slc22a30", "Slc7a13", "Bcat1", "Slc2a1"),
  Profibrotic_PT = c("Pdgfb", "Cd74", "Tnfrsf12a", "Cxcl1", "Cxcl10", "Cxcl16"),
  Transient_mix  = c("Wfdc15b"),
  Proliferating  = c("Mki67"),
  Immune         = c("Cd52")
)

# remove the old PT_score from script 2
old_score_cols <- grep("_score1$", colnames(ckd_pt_uuo@meta.data), value = TRUE)
ckd_pt_uuo@meta.data[, old_score_cols] <- NULL

# for each marker gene set
for (name in names(pt_marker_sets)) {
  genes <- pt_marker_sets[[name]]
  genes <- genes[genes %in% rownames(ckd_pt_uuo)]
  
  # check how strongly each cell express this set of genes.
  ckd_pt_uuo <- AddModuleScore(
    ckd_pt_uuo,
    features = list(genes),
    name = paste0(name, "_score")
  )
}

score_cols <- grep("_score1$", colnames(ckd_pt_uuo@meta.data), value = TRUE)

# average the score across all cells in each cluster
# this gives one mean PT score per cluster
cluster_scores <- aggregate(
  ckd_pt_uuo@meta.data[, score_cols],
  by = list(cluster = ckd_pt_uuo$seurat_clusters),
  FUN = mean
)

# count number of cells in each cluster
cluster_sizes <- as.data.frame(table(ckd_pt_uuo$seurat_clusters))
colnames(cluster_sizes) <- c("cluster", "n_cells")
cluster_sizes$cluster <- as.character(cluster_sizes$cluster)

# make score table easier to work with
cluster_scores$cluster <- as.character(cluster_scores$cluster)
score_only <- cluster_scores[, -1, drop = FALSE]

# make a quick summary table
compact_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(cluster_scores)), function(i) {
    cluster_id <- cluster_scores$cluster[i]
    scores <- as.numeric(score_only[i, ])
    names(scores) <- colnames(score_only)
    
    ord <- order(scores, decreasing = TRUE)
    
    data.frame(
      cluster = cluster_id,
      top_label = sub("_score1$", "", names(scores)[ord[1]]),
      top_score = scores[ord[1]],
      second_label = sub("_score1$", "", names(scores)[ord[2]]),
      second_score = scores[ord[2]],
      profibrotic_score = cluster_scores$Profibrotic_PT_score1[i],
      row.names = NULL
    )
  })
)

compact_summary <- merge(compact_summary, cluster_sizes, by = "cluster")
compact_summary <- compact_summary[, c(
  "cluster", "n_cells",
  "top_label", "top_score",
  "second_label", "second_score",
  "profibrotic_score"
)]

compact_summary <- compact_summary[order(as.numeric(compact_summary$cluster)), ]
compact_summary
save_table(compact_summary, "03_PT_UUO_label_scores.tsv")

# label each cluster w/ top scoring labels

# create cluster to label mapping and assign labels to each cell.
# create cluster → label mapping
cluster_to_label <- setNames(
  compact_summary$top_label,
  compact_summary$cluster
)

# make sure clusters are character
clusters_char <- as.character(ckd_pt_uuo$seurat_clusters)

# map cluster → label (no names attached)
labels_vec <- cluster_to_label[clusters_char]
labels_vec <- as.vector(labels_vec)

# assign to metadata
ckd_pt_uuo$pt_ident_pred <- labels_vec

# plot labeled UMAP
p_umap_pt_labels <- DimPlot(
  ckd_pt_uuo,
  reduction = "umap",
  group.by = "pt_ident_pred",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) + NoLegend()

p_umap_pt_labels

save_plot(
  p_umap_pt_labels,
  "03_PT_UUO_umap_pt_ident_pred.png"
)

saveRDS(ckd_pt_uuo,
        file = file.path(out_obj_dir, "03_gse182256_labeled_PT.rds"))