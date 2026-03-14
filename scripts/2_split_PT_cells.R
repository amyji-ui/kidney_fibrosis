# Qinglin Kong
# 3/13/26

# Load processed Seurat object
# Find cluster marker genes
# Identify proximal tubule (PT) clusters using canonical PT markers
# Split the dataset into PT and non-PT subsets
# Save PT and non-PT Seurat objects for downstream analysis

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

ckd <- readRDS("results/seurat_objects/01_gse182256_processed_umap.rds")

# -----------------
# Find differentially expressed features (cluster markers)
# -----------------

# find all differentially expressed markers
ckd_markers <- FindAllMarkers(ckd, only.pos = TRUE, logfc.threshold = 0.25)
top10_markers <- ckd_markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 10)

# save the markers to file
dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)

write.table(
  ckd_markers,
  file = "results/markers/02_all_cluster_markers.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top10_markers,
  file = "results/markers/02_top10_markers_per_cluster.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------
# Use PT markers to find clusters that represent PT cells
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6683720/
# Table 1 says the gold standard markers for All proximal tubule are
# Slc34a1 (Na-Pi2) and Lrp2 (megalin)
# -----------------

# set cluster identities so Seurat functions work on each cluster
Idents(ckd) <- "seurat_clusters"
cluster_ids <- levels(Idents(ckd))

# for each cluster, test if Lrp2 and Slc34a1 are enriched relative to the rest
# of the dataset.
pt_lrp2 <- lapply(
  cluster_ids,
  function(cl) {
    FindMarkers(ckd, ident.1 = cl, features = "Lrp2")
  }
)
pt_slc34a1 <- lapply(
  cluster_ids,
  function(cl) {
    FindMarkers(ckd, ident.1 = cl, features = "Slc34a1")
  }
)

names(pt_lrp2) <- cluster_ids
names(pt_slc34a1) <- cluster_ids

# build a summary table of marker statistics for each cluster
pt_stats <- data.frame(
  cluster = cluster_ids,

  Lrp2_padj = sapply(
    pt_lrp2,
    function(marker_result) {
      marker_result["Lrp2", "p_val_adj"]
    }
  ),

  Slc34a1_padj = sapply(
    pt_slc34a1,
    function(marker_result) {
      marker_result["Slc34a1", "p_val_adj"]
    }
  )
)
pt_stats

# compute a PT gene module score
# by combining the two PT marker gene into a simple module score to measure
# overall PT-like transcriptional activity per cell,
pt_genes <- c("Lrp2", "Slc34a1")
pt_genes <- pt_genes[pt_genes %in% rownames(ckd)]
ckd <- AddModuleScore(ckd, features = list(pt_genes), name = "PT_score")

# summarize the average PT score per cluster
pt_score_by_cluster <- aggregate(ckd$PT_score1,
                                 by = list(cluster = ckd$seurat_clusters),
                                 FUN = mean)
colnames(pt_score_by_cluster)[2] <- "mean_PT_score"

pt_expr <- FetchData(ckd, vars = c("seurat_clusters", "Lrp2", "Slc34a1"))

# calculate fraction of cells expressing PT markers
pt_pct_summary <- do.call(
  rbind,
  lapply(split(pt_expr, pt_expr$seurat_clusters), function(df) {
    data.frame(
      cluster = unique(as.character(df$seurat_clusters)),
      pct_Lrp2 = mean(df$Lrp2 > 0),
      pct_Slc34a1 = mean(df$Slc34a1 > 0)
    )
  })
)

pt_summary <- merge(pt_stats, pt_score_by_cluster, by = "cluster")
pt_summary <- merge(pt_summary, pt_pct_summary, by = "cluster")

# sort clusters by PT score for easier inspection
pt_summary <- pt_summary[order(-pt_summary$mean_PT_score), ]
pt_summary

# save the table
write.table(
  pt_summary,
  file = file.path(out_obj_dir, "02_pt_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# define PT clusters:
# 1) both canonical PT markers significantly enriched
# 2) both markers have positive effect size
# 3) a substantial fraction of cells express both markers
pt_clusters <- pt_summary$cluster[
  !is.na(pt_summary$Lrp2_padj) &
    !is.na(pt_summary$Slc34a1_padj) &
    pt_summary$Lrp2_padj < 0.05 &
    pt_summary$Slc34a1_padj < 0.05 &
    pt_summary$mean_PT_score > 0 &
    (pt_summary$pct_Lrp2 > 0.20 | pt_summary$pct_Slc34a1 > 0.20)
]
pt_clusters

# -----------------
# Subset PT clusters
# -----------------
ckd_pt <- subset(ckd, idents = pt_clusters)
ckd_nonpt <- subset(ckd, idents = pt_clusters, invert = TRUE)
ncol(ckd)
ncol(ckd_pt)
ncol(ckd_nonpt)

table(ckd_pt$seurat_clusters)
table(ckd_nonpt$seurat_clusters)

saveRDS(ckd_pt,
        file = file.path(out_obj_dir, "02_gse182256_PT_subset.rds"))
saveRDS(ckd_nonpt,
        file = file.path(out_obj_dir, "02_gse182256_nonPT_subset.rds"))
