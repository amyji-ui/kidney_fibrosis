library(Matrix)
library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)

mat_path  <- "data/GSE182256_Export_counts.rds"
meta_path <- "data/GSE182256_Export_Metadata.txt"
umap_path <- "data/GSE182256_EXPORT_umap.txt"

mat <- readRDS(mat_path)
obj <- CreateSeuratObject(counts = mat)

meta <- fread(meta_path, sep = "\t",
              header = TRUE,
              fill = TRUE,
              quote = "",
              data.table = FALSE)
names(meta)[1] <- "cell"

umap <- fread(umap_path,
              sep = "\t",
              header = TRUE,
              fill = TRUE,
              quote = "",
              data.table = FALSE)
names(umap)[1] <- "cell"

# align metadata to Seurat cells
rownames(meta) <- meta$cell

# keep only overlapping cells
common <- intersect(colnames(obj), rownames(meta))
obj <- subset(obj, cells = common)
obj <- AddMetaData(obj, metadata = meta)

# confirm Cluster now exists
head(obj@meta.data$Cluster)

#============================== Figure 1b ==============================
df <- merge(umap, meta[, c("cell", "Cluster")], by = "cell")

# numeric labels for each cell type
ct_levels <- sort(unique(as.character(df$Cluster)))
df$ClusterNum <- match(as.character(df$Cluster), ct_levels)

# cluster centers (median position) for label placement
centers <- aggregate(df[, c("UMAP_1", "UMAP_2")],
                     by = list(ClusterNum = df$ClusterNum),
                     FUN = median)
# plot
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(size = 0.15, alpha = 0.8) +
  geom_text(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = ClusterNum),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_classic()

ggsave("output/figure1/fig1b_umap.png", plot = p, width = 7, height = 5, dpi = 300)
mapping_df <- data.frame(Number = seq_along(ct_levels), CellType = ct_levels)
write.csv(mapping_df, "output/figure1/cluster_number_mapping.csv", row.names = FALSE)

# ToDo: the above figure1b was generated using the given UMAP data.
# We should attempt to recreate this from the count matrix.

#========================== figure 1c ===============================
# find the top genes.
obj$Condition <- ifelse(grepl("^UUO", obj$orig.ident), "UUO", "Sham")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- obj$Cluster

# Not sure if the sparse matrix is normalized or not. Normalize just in case. 
obj <- NormalizeData(obj)

# Finds markers (differentially expressed genes) for each of the identity classes in a dataset
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# pick top 1 genes per cluster (adjust n=1..3)
# Need a while to compute! 
top <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1) %>%
  ungroup()

#store the result
genes <- unique(top$gene)
genes_df <- data.frame(gene = genes)
write.csv(genes_df, "output/figure1/fig1c_genes.csv", row.names = FALSE)

# Generate plot
obj$Condition <- ifelse(grepl("^UUO", obj$orig.ident), "UUO", "Sham")
obj$Condition <- factor(obj$Condition, levels = c("UUO", "Sham"))

# Make Seurat DotPlot object (we'll steal its computed data)
dp <- DotPlot(
  obj,
  features = genes,
  group.by = "Cluster",
  split.by = "Condition"
)

df <- dp$data

# df$id looks like "A-IC_UUO" etc. Split it into CellType + Condition
df <- df %>%
  mutate(
    CellType = sub("_(UUO|Sham)$", "", id),
    Condition = sub("^.*_(UUO|Sham)$", "\\1", id),
    Condition = factor(Condition, levels = c("UUO", "Sham"))
  )

#use expression strength as transparency (keeps red/blue color by condition)
# avg.exp.scaled is typically in the DotPlot data; if not, use avg.exp
if ("avg.exp.scaled" %in% names(df)) {
  rng <- range(df$avg.exp.scaled, na.rm = TRUE)
  df$alpha_expr <- (df$avg.exp.scaled - rng[1]) / (rng[2] - rng[1] + 1e-9)
  df$alpha_expr <- 0.2 + 0.8 * df$alpha_expr  # keep visible
} else {
  df$alpha_expr <- 0.8
}

p <- ggplot(df, aes(x = features.plot, y = CellType)) +
  geom_point(
    aes(size = pct.exp, color = Condition, alpha = alpha_expr),
    position = position_dodge(width = 0.6)
  ) +
  scale_color_manual(values = c(UUO = "red", Sham = "blue")) +
  scale_alpha_identity() +
  labs(x = NULL, y = NULL, size = "% expressing") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/figure1/fig1c_1.png", p, width = 14, height = 6, dpi = 300)