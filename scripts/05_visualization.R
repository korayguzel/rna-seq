#!/usr/bin/env Rscript
# Publication-quality visualizations for mouse brain aging RNA-seq
# GSE225576 - Brain tissue only

message("Loading packages...")
library(tidyverse)
library(pheatmap)

# Check if EnhancedVolcano is available
ev_available <- requireNamespace("EnhancedVolcano", quietly = TRUE)
if (ev_available) {
  library(EnhancedVolcano)
  message("  Using EnhancedVolcano for volcano plot")
} else {
  message("  EnhancedVolcano not available, using ggplot2 for volcano plot")
}

# Create figures directory if needed
if (!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}

# Load data
message("\nLoading data...")
de_results <- read.csv("results/tables/DE_30m_vs_6m.csv")
vst_matrix <- readRDS("data/processed/brain_vst.rds")
metadata <- readRDS("data/processed/brain_metadata.rds")

message("  DE results: ", nrow(de_results), " genes")
message("  VST matrix: ", nrow(vst_matrix), " genes x ", ncol(vst_matrix), " samples")

# ==============================================================================
# 1. VOLCANO PLOT (30m vs 6m)
# ==============================================================================
message("\nGenerating volcano plot...")

# Define significance thresholds
lfc_cutoff <- 1
padj_cutoff <- 0.05

# Add significance column
de_results$significance <- "NS"
de_results$significance[de_results$log2FoldChange > lfc_cutoff &
                         de_results$padj < padj_cutoff] <- "Up"
de_results$significance[de_results$log2FoldChange < -lfc_cutoff &
                         de_results$padj < padj_cutoff] <- "Down"
de_results$significance <- factor(de_results$significance, levels = c("Up", "Down", "NS"))

# Get top 10 genes by padj for labeling
top_genes <- head(de_results[order(de_results$padj), ], 10)

if (ev_available) {
  # Use EnhancedVolcano
  p_volcano <- EnhancedVolcano(
    de_results,
    lab = de_results$gene_symbol,
    x = "log2FoldChange",
    y = "pvalue",  # Uses raw p-value, will convert to -log10 internally
    title = "Differential Expression: 30m vs 6m",
    subtitle = "Brain aging RNA-seq",
    caption = paste0("Total genes: ", nrow(de_results)),
    pCutoff = padj_cutoff,
    FCcutoff = lfc_cutoff,
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.8,
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    legendPosition = "right",
    legendLabSize = 10,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    boxedLabels = FALSE,
    max.overlaps = 15
  )

  ggsave("results/figures/volcano_30m_vs_6m.png", p_volcano,
         dpi = 300, width = 10, height = 8, units = "in")
} else {
  # Use ggplot2
  p_volcano <- ggplot(de_results, aes(x = log2FoldChange, y = -log10(pvalue),
                                       color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey60"),
                       labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not significant")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text(data = top_genes, aes(label = gene_symbol),
              size = 3, vjust = -1, hjust = 0.5, color = "black") +
    labs(x = expression(Log[2]~Fold~Change),
         y = expression(-Log[10]~P-value),
         title = "Differential Expression: 30m vs 6m",
         subtitle = "Brain aging RNA-seq",
         color = "Significance") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")

  ggsave("results/figures/volcano_30m_vs_6m.png", p_volcano,
         dpi = 300, width = 10, height = 8, units = "in", bg = "white")
}

message("  Saved: results/figures/volcano_30m_vs_6m.png")

# ==============================================================================
# 2. HEATMAP (Top 50 DE genes)
# ==============================================================================
message("\nGenerating heatmap...")

# Get top 50 genes by padj
top50_genes <- head(de_results[order(de_results$padj), ], 50)$gene_id

# Extract VST data for top 50 genes
heatmap_data <- vst_matrix[rownames(vst_matrix) %in% top50_genes, ]
heatmap_data <- heatmap_data[match(top50_genes, rownames(heatmap_data)), ]

# Ensure column order matches metadata
heatmap_data <- heatmap_data[, metadata$sample_id]

# Create annotation
annotation_col <- data.frame(
  Age = metadata$age_group,
  row.names = metadata$sample_id
)

annotation_colors <- list(
  Age = c("6m" = "#1f77b4", "15m" = "#ff7f0e", "24m" = "#2ca02c", "30m" = "#d62728")
)

png("results/figures/heatmap_top50.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Top 50 Differentially Expressed Genes",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

message("  Saved: results/figures/heatmap_top50.png")

# ==============================================================================
# 2b. HEATMAP WITH GENE CLUSTER ANNOTATION (Prompt 5.1)
# ==============================================================================
message("\nGenerating annotated heatmap with gene clusters...")

# Define gene categories
gene_categories <- list(
  "Immune/Inflammation" = c("Cd74", "Lyz2", "Clec7a", "C3", "C4b", "Ctss", "Stat1",
                            "Gpnmb", "Lgals3", "Ifit3b", "Oasl2", "H2-Q4", "H2-Q6",
                            "H2-K1", "H2-D1", "B2m"),
  "Synaptic/Neuronal" = c("Syp", "Snap25", "Camk2a", "Ank1"),
  "Astrocyte" = c("Gfap", "Vim")
)

# Assign categories to top 50 genes
gene_cluster <- data.frame(
  gene_id = top50_genes,
  cluster = "Others",
  stringsAsFactors = FALSE
)

for (category in names(gene_categories)) {
  genes_in_category <- gene_categories[[category]]
  gene_cluster$cluster[gene_cluster$gene_id %in% genes_in_category] <- category
}

# Create row annotation for pheatmap
annotation_row <- data.frame(
  Cluster = factor(gene_cluster$cluster,
                   levels = c("Immune/Inflammation", "Synaptic/Neuronal", "Astrocyte", "Others")),
  row.names = gene_cluster$gene_id
)

# Define annotation colors
annotation_colors_row <- list(
  Cluster = c("Immune/Inflammation" = "#e41a1c",
              "Synaptic/Neuronal" = "#377eb8",
              "Astrocyte" = "#4daf4a",
              "Others" = "#999999"),
  Age = c("6m" = "#1f77b4", "15m" = "#ff7f0e", "24m" = "#2ca02c", "30m" = "#d62728")
)

png("results/figures/heatmap_top50_annotated.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors_row,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Top 50 DE Genes with Functional Clusters",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

message("  Saved: results/figures/heatmap_top50_annotated.png")

# Print cluster summary
cluster_summary <- table(gene_cluster$cluster)
message("  Gene cluster distribution:")
for (cl in names(cluster_summary)) {
  message("    ", cl, ": ", cluster_summary[cl], " genes")
}

# ==============================================================================
# 2c. POLISHED VOLCANO PLOT WITH GGREPEL (Prompt 5.2)
# ==============================================================================
message("\nGenerating polished volcano plot with ggrepel...")

# Check/install ggrepel
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("  Installing ggrepel...")
  install.packages("ggrepel", repos = "https://cloud.r-project.org/")
}
library(ggrepel)

# Get top 15 genes by padj for labeling
top15_genes <- head(de_results[order(de_results$padj), ], 15)

# Genes to force label (must be visible)
force_genes <- c("Gfap", "C4b", "Clec7a", "C3", "Gm2446")
force_data <- de_results[de_results$gene_symbol %in% force_genes, ]

# Combine: top 15 + forced genes (without duplicates)
label_genes <- unique(c(top15_genes$gene_symbol, force_data$gene_symbol))
label_data <- de_results[de_results$gene_symbol %in% label_genes, ]

message("  Labeling ", nrow(label_data), " genes (top 15 + forced)")

# Create polished volcano plot
p_volcano_polished <- ggplot(de_results, aes(x = log2FoldChange, y = -log10(pvalue),
                                              color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey60"),
                     labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not significant")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = label_data,
                  aes(label = gene_symbol),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.2,
                  color = "black",
                  force = 10) +
  labs(x = expression(Log[2]~Fold~Change),
       y = expression(-Log[10]~P-value),
       title = "Differential Expression: 30m vs 6m",
       subtitle = "Brain aging RNA-seq - Polished with ggrepel",
       color = "Significance") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

ggsave("results/figures/volcano_30m_vs_6m_polished.png", p_volcano_polished,
       dpi = 300, width = 10, height = 8, units = "in", bg = "white")

message("  Saved: results/figures/volcano_30m_vs_6m_polished.png")

# ==============================================================================
# 3. MARKER GENE BOXPLOTS
# ==============================================================================
message("\nGenerating marker gene boxplots...")

# Define marker genes
up_markers <- c("Gfap", "C4b", "Ctss")
down_markers <- c("Syp", "Snap25", "Camk2a")
all_markers <- c(up_markers, down_markers)

# Check which markers are in our data
available_markers <- all_markers[all_markers %in% rownames(vst_matrix)]
missing_markers <- all_markers[!all_markers %in% rownames(vst_matrix)]

if (length(missing_markers) > 0) {
  message("  Warning: Missing markers: ", paste(missing_markers, collapse = ", "))
}

if (length(available_markers) == 0) {
  message("  Error: No marker genes found in dataset!")
} else {
  # Prepare data for plotting
  plot_data <- data.frame()

  for (gene in available_markers) {
    gene_data <- data.frame(
      gene = gene,
      sample_id = colnames(vst_matrix),
      expression = vst_matrix[gene, ],
      age_group = metadata$age_group[match(colnames(vst_matrix), metadata$sample_id)]
    )
    plot_data <- rbind(plot_data, gene_data)
  }

  # Set gene order and category
  plot_data$category <- ifelse(plot_data$gene %in% up_markers, "Up (Neuroinflammation)", "Down (Synaptic)")
  plot_data$gene <- factor(plot_data$gene, levels = available_markers)
  plot_data$age_group <- factor(plot_data$age_group, levels = c("6m", "15m", "24m", "30m"))

  # Create plot
  p_markers <- ggplot(plot_data, aes(x = age_group, y = expression, fill = age_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("6m" = "#1f77b4", "15m" = "#ff7f0e",
                                  "24m" = "#2ca02c", "30m" = "#d62728")) +
    labs(x = "Age Group", y = "VST Normalized Expression",
         title = "Aging Marker Genes Across Age Groups") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 0),
          legend.position = "none")

  ggsave("results/figures/marker_genes.png", p_markers,
         dpi = 300, width = 12, height = 6, units = "in", bg = "white")

  message("  Saved: results/figures/marker_genes.png")
  message("  Plotted markers: ", paste(available_markers, collapse = ", "))

  # ==============================================================================
  # 3b. MARKER GENES WITH STATISTICAL SIGNIFICANCE (Prompt 5.3)
  # ==============================================================================
  message("\nGenerating marker genes with significance testing...")

  # Check/install ggsignif
  if (!requireNamespace("ggsignif", quietly = TRUE)) {
    message("  Installing ggsignif...")
    install.packages("ggsignif", repos = "https://cloud.r-project.org/")
  }
  library(ggsignif)

  # Function to get significance stars
  get_stars <- function(p) {
    if (is.na(p) || p > 0.05) return("ns")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    return("*")
  }

  # Calculate Wilcoxon test for each gene (6m vs 30m only)
  message("\n  Wilcoxon test results (6m vs 30m):")

  sig_data <- data.frame()

  for (gene in available_markers) {
    gene_data <- plot_data[plot_data$gene == gene & plot_data$age_group %in% c("6m", "30m"), ]

    # Wilcoxon rank-sum test
    test_result <- wilcox.test(expression ~ age_group, data = gene_data)
    pval <- test_result$p.value
    stars <- get_stars(pval)

    message(sprintf("    %s: p = %.4f (%s)", gene, pval, stars))

    # Get y position for bracket
    y_max <- max(plot_data$expression[plot_data$gene == gene], na.rm = TRUE)
    y_pos <- y_max + 0.5

    sig_data <- rbind(sig_data, data.frame(
      gene = gene,
      p_value = pval,
      stars = stars,
      y_position = y_pos,
      x1 = "6m",
      x2 = "30m"
    ))
  }

  # Prepare bracket data for geom_segment (proper facet-aware annotations)
  bracket_data <- data.frame()
  for (i in 1:nrow(sig_data)) {
    if (sig_data$stars[i] == "ns") next

    bracket_data <- rbind(bracket_data, data.frame(
      gene = sig_data$gene[i],
      x = c(1, 1, 4),  # horizontal line, left tick, right tick x positions
      xend = c(4, 1, 4),
      y = c(sig_data$y_position[i], sig_data$y_position[i], sig_data$y_position[i]),
      yend = c(sig_data$y_position[i], sig_data$y_position[i] - 0.1, sig_data$y_position[i] - 0.1)
    ))
  }

  # Text labels for stars
  text_data <- sig_data[sig_data$stars != "ns", ]

  # Create plot with significance brackets using geom_segment (facet-aware)
  p_markers_sig <- ggplot(plot_data, aes(x = age_group, y = expression, fill = age_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("6m" = "#1f77b4", "15m" = "#ff7f0e",
                                  "24m" = "#2ca02c", "30m" = "#d62728")) +
    labs(x = "Age Group", y = "VST Normalized Expression",
         title = "Aging Marker Genes - Statistical Significance",
         subtitle = "Wilcoxon test: 6m vs 30m") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(angle = 0),
          legend.position = "none")

  # Add brackets only if there are significant genes
  if (nrow(bracket_data) > 0) {
    p_markers_sig <- p_markers_sig +
      geom_segment(data = bracket_data, aes(x = x, xend = xend, y = y, yend = yend),
                   color = "black", linewidth = 0.5, inherit.aes = FALSE)
  }

  # Add star labels
  if (nrow(text_data) > 0) {
    p_markers_sig <- p_markers_sig +
      geom_text(data = text_data, aes(x = 2.5, y = y_position + 0.2, label = stars),
                size = 4, fontface = "bold", color = "black", inherit.aes = FALSE)
  }

  ggsave("results/figures/marker_genes_significance.png", p_markers_sig,
         dpi = 300, width = 12, height = 6, units = "in", bg = "white")

  message("\n  Saved: results/figures/marker_genes_significance.png")
}

# ==============================================================================
# 4. DE SUMMARY BARPLOT
# ==============================================================================
message("\nGenerating DE summary barplot...")

# Load all DE results
comparisons <- c("15m_vs_6m", "24m_vs_15m", "30m_vs_24m", "30m_vs_6m")
summary_data <- data.frame()

for (comp in comparisons) {
  de_file <- paste0("results/tables/DE_", comp, ".csv")
  if (file.exists(de_file)) {
    de_res <- read.csv(de_file)

    # Count significant genes
    sig_up <- sum(!is.na(de_res$padj) & de_res$padj < 0.05 & de_res$log2FoldChange > 1, na.rm = TRUE)
    sig_down <- sum(!is.na(de_res$padj) & de_res$padj < 0.05 & de_res$log2FoldChange < -1, na.rm = TRUE)

    summary_data <- rbind(summary_data, data.frame(
      comparison = comp,
      direction = "Up",
      count = sig_up
    ))
    summary_data <- rbind(summary_data, data.frame(
      comparison = comp,
      direction = "Down",
      count = sig_down
    ))
  }
}

# Format comparison labels
summary_data$comparison <- factor(summary_data$comparison,
                                   levels = comparisons,
                                   labels = c("15m vs 6m", "24m vs 15m",
                                             "30m vs 24m", "30m vs 6m"))
summary_data$direction <- factor(summary_data$direction, levels = c("Up", "Down"))

# Create barplot
p_summary <- ggplot(summary_data, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"),
                    labels = c("Up" = "Upregulated", "Down" = "Downregulated")) +
  labs(x = "Comparison", y = "Number of DEGs",
       title = "Differentially Expressed Genes Summary",
       subtitle = "padj < 0.05 & |log2FC| > 1",
       fill = "Direction") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave("results/figures/DE_summary.png", p_summary,
       dpi = 300, width = 8, height = 6, units = "in", bg = "white")

message("  Saved: results/figures/DE_summary.png")

# ==============================================================================
# Complete
# ==============================================================================
message("\n========================================")
message("All visualizations generated successfully!")
message("========================================")
message("\nGenerated files:")
message("  - results/figures/volcano_30m_vs_6m.png")
message("  - results/figures/volcano_30m_vs_6m_polished.png (Prompt 5.2)")
message("  - results/figures/heatmap_top50.png")
message("  - results/figures/heatmap_top50_annotated.png (Prompt 5.1)")
message("  - results/figures/marker_genes.png")
message("  - results/figures/marker_genes_significance.png (Prompt 5.3)")
message("  - results/figures/DE_summary.png")
