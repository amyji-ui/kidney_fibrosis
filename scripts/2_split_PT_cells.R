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
  # for cluster cl, it tests expression in cluster cl vs all other cells
  # returns avg_log2FC which is how much higher the exp is
  # also returns p_val_adj which is the statistical significance
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

# name each list element by cluster ID so results are easier to inspect
names(pt_lrp2) <- cluster_ids
names(pt_slc34a1) <- cluster_ids

# make a small summary table that keeps only p_val_adj
# for the two canonical PT markers in each cluster
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

# define the PT marker gene set used for a simple PT module score
pt_genes <- c("Lrp2", "Slc34a1")
pt_genes <- pt_genes[pt_genes %in% rownames(ckd)]

# add a PT module score to each cell
# higher score means more PT-like expression pattern
ckd <- AddModuleScore(ckd, features = list(pt_genes), name = "PT_score")

# average the PT module score across all cells in each cluster
# this gives one mean PT score per cluster
pt_score_by_cluster <- aggregate(ckd$PT_score1,
                                 by = list(cluster = ckd$seurat_clusters),
                                 FUN = mean)
colnames(pt_score_by_cluster)[2] <- "mean_PT_score"

# pull out per-cell expression values for the two PT markers
# together with cluster identity
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
# 2) a substantial fraction of cells express both markers
# 3) at least 20% cells express these markers in the cluster
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
