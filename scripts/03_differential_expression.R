#!/usr/bin/env Rscript
# Differential expression analysis with DESeq2
# GSE225576 - Mouse brain aging RNA-seq

message("Loading packages...")
library(DESeq2)
library(tidyverse)

# Create results directory if needed
if (!dir.exists("results/tables")) {
  dir.create("results/tables", recursive = TRUE)
}

message("Loading data...")
counts_filtered <- readRDS("data/processed/brain_counts_filtered.rds")
metadata <- readRDS("data/processed/brain_metadata.rds")

message("Samples: ", ncol(counts_filtered))
message("Genes: ", nrow(counts_filtered))

# Ensure counts are integers
if (!all(counts_filtered == floor(counts_filtered), na.rm = TRUE)) {
  message("Rounding counts to integers...")
  counts_filtered <- round(counts_filtered)
}

# 1. Create DESeqDataSet
message("\nCreating DESeqDataSet...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ age_group
)

# 2. Run DESeq
message("Running DESeq...")
dds <- DESeq(dds)

# Get results names to see available contrasts
message("\nAvailable contrasts:")
resultsNames(dds) %>% print()

# 3. Define comparisons
comparisons <- list(
  list(name = "15m_vs_6m", contrast = c("age_group", "15m", "6m")),
  list(name = "24m_vs_15m", contrast = c("age_group", "24m", "15m")),
  list(name = "30m_vs_24m", contrast = c("age_group", "30m", "24m")),
  list(name = "30m_vs_6m", contrast = c("age_group", "30m", "6m"))
)

# 4. Process each comparison
for (comp in comparisons) {
  comp_name <- comp$name
  contrast <- comp$contrast

  message("\n", paste(rep("=", 50), collapse = ""))
  message("Comparison: ", contrast[2], " vs ", contrast[3])
  message(paste(rep("=", 50), collapse = ""))

  # Get results and shrink LFC
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res_shrunk <- lfcShrink(dds, contrast = contrast, type = "ashr", res = res)

  # Convert to data frame
  res_df <- as.data.frame(res_shrunk)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_symbol <- res_df$gene_id  # Gene IDs are already symbols
  rownames(res_df) <- NULL

  # Add stat column if missing (lfcShrink doesn't return it)
  if (!"stat" %in% names(res_df)) {
    res_df$stat <- NA
  }

  # Reorder columns
  res_df <- res_df[, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange",
                       "lfcSE", "stat", "pvalue", "padj")]

  # Filter: padj < 0.05 AND |log2FoldChange| > 1
  sig_genes <- res_df[!is.na(res_df$padj) &
                        res_df$padj < 0.05 &
                        abs(res_df$log2FoldChange) > 1, ]

  up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
  down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

  message("Upregulated: ", nrow(up_genes), " genes")
  message("Downregulated: ", nrow(down_genes), " genes")

  # Top 5 upregulated
  if (nrow(up_genes) > 0) {
    top_up <- head(up_genes[order(up_genes$log2FoldChange, decreasing = TRUE), ], 5)
    message("Top 5 up: ", paste(top_up$gene_symbol, collapse = ", "))
  }

  # Top 5 downregulated
  if (nrow(down_genes) > 0) {
    top_down <- head(down_genes[order(down_genes$log2FoldChange), ], 5)
    message("Top 5 down: ", paste(top_down$gene_symbol, collapse = ", "))
  }

  # Save full results table
  output_file <- paste0("results/tables/DE_", comp_name, ".csv")
  write.csv(res_df, output_file, row.names = FALSE)
  message("Saved: ", output_file)
}

message("\n", paste(rep("=", 50), collapse = ""))
message("Differential expression analysis complete.")
message(paste(rep("=", 50), collapse = ""))

message("\nSession info:")
sessionInfo()
