#!/usr/bin/env Rscript
# QC and filtering for mouse brain aging RNA-seq
# GSE225576 - Brain tissue only

message("Loading packages...")
library(tidyverse)
library(DESeq2)
library(pheatmap)

# Create results directory if needed
if (!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}

message("Loading data...")
counts <- readRDS("data/processed/brain_counts.rds")
metadata <- readRDS("data/processed/brain_metadata.rds")

message("Samples: ", ncol(counts))
message("Genes: ", nrow(counts))

# 1. Library Size Plot
message("\nCreating library size plot...")
library_sizes <- colSums(counts)
lib_df <- data.frame(
  sample_id = names(library_sizes),
  total_counts = library_sizes,
  age_group = metadata$age_group[match(names(library_sizes), metadata$sample_id)]
)

p_lib <- ggplot(lib_df, aes(x = sample_id, y = total_counts / 1e6, fill = age_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("6m" = "#1f77b4", "15m" = "#ff7f0e",
                               "24m" = "#2ca02c", "30m" = "#d62728")) +
  labs(x = "Sample", y = "Total Counts (Millions)",
       title = "Library Size by Sample", fill = "Age Group") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/qc_library_size.png", p_lib,
       dpi = 300, width = 8, height = 6, units = "in")
message("  Saved: results/figures/qc_library_size.png")

# 2. Gene Filtering
message("\nFiltering genes by CPM...")

# Calculate CPM manually: (counts / library_size) * 1e6
lib_sizes <- colSums(counts)
cpm <- sweep(counts, 2, lib_sizes, FUN = "/") * 1e6
# Keep genes with CPM > 1 in at least 3 samples
keep <- rowSums(cpm > 1) >= 3

message("  Before filtering: ", nrow(counts), " genes")
counts_filtered <- counts[keep, ]
message("  After filtering: ", nrow(counts_filtered), " genes")

# 3. Normalization with VST
message("\nRunning VST normalization...")

# DESeq2 requires integer counts
if (!all(counts_filtered == floor(counts_filtered), na.rm = TRUE)) {
  message("  Rounding counts to integers for DESeq2...")
  counts_filtered <- round(counts_filtered)
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ age_group
)

vst_data <- vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)

# 4. PCA Plot
message("\nCreating PCA plot...")

# Select top 500 most variable genes
row_vars <- apply(vst_matrix, 1, var)
top_genes <- names(sort(row_vars, decreasing = TRUE))[1:500]
vst_top <- vst_matrix[top_genes, ]

# Run PCA
pca <- prcomp(t(vst_top))
pca_var <- summary(pca)$importance
pc1_var <- round(pca_var[2, 1] * 100, 1)
pc2_var <- round(pca_var[2, 2] * 100, 1)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  sample_id = rownames(pca$x),
  age_group = metadata$age_group[match(rownames(pca$x), metadata$sample_id)]
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = age_group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("6m" = "#1f77b4", "15m" = "#ff7f0e",
                                "24m" = "#2ca02c", "30m" = "#d62728")) +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       title = "PCA of Brain Samples (Top 500 Variable Genes)",
       color = "Age Group") +
  theme_bw(base_size = 12)

ggsave("results/figures/qc_pca.png", p_pca,
       dpi = 300, width = 8, height = 8, units = "in")
message("  Saved: results/figures/qc_pca.png")

# 5. Correlation Heatmap
message("\nCreating correlation heatmap...")

# Pearson correlation between samples
cor_matrix <- cor(vst_matrix, method = "pearson")

# Annotation for heatmap
anno_row <- data.frame(
  age_group = metadata$age_group,
  row.names = metadata$sample_id
)

anno_colors <- list(age_group = c("6m" = "#1f77b4", "15m" = "#ff7f0e",
                                   "24m" = "#2ca02c", "30m" = "#d62728"))

png("results/figures/qc_correlation.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(cor_matrix,
         annotation_row = anno_row,
         annotation_col = anno_row,
         annotation_colors = anno_colors,
         main = "Sample Correlation Heatmap",
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize = 8)
dev.off()
message("  Saved: results/figures/qc_correlation.png")

# Save outputs
message("\nSaving processed data...")
saveRDS(counts_filtered, "data/processed/brain_counts_filtered.rds")
saveRDS(vst_matrix, "data/processed/brain_vst.rds")

message("  Saved: data/processed/brain_counts_filtered.rds")
message("  Saved: data/processed/brain_vst.rds")

message("\nQC complete.")
