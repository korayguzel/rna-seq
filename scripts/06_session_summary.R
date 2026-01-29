#!/usr/bin/env Rscript
# Session info and project summary
# GSE225576 - Mouse brain aging RNA-seq

message("Loading packages...")
library(tidyverse)

# ==============================================================================
# 1. SAVE SESSION INFO
# ==============================================================================
message("\nStep 1: Saving session info...")

# Capture session info
si <- capture.output(sessionInfo())

# Write to file
writeLines(si, "results/session_info.txt")
message("  Session info saved to results/session_info.txt")

# ==============================================================================
# 2. CREATE ANALYSIS SUMMARY TABLE
# ==============================================================================
message("\nStep 2: Creating analysis summary table...")

# Load necessary data
metadata <- readRDS("data/processed/brain_metadata.rds")
counts_raw <- readRDS("data/processed/brain_counts.rds")
counts_filtered <- readRDS("data/processed/brain_counts_filtered.rds")
de_results <- read.csv("results/tables/DE_30m_vs_6m.csv")

# Calculate sample counts per group
samples_per_group <- table(metadata$age_group)

# Get DE counts for all comparisons
comparisons <- c("15m_vs_6m", "24m_vs_15m", "30m_vs_24m", "30m_vs_6m")
de_counts <- data.frame()

for (comp in comparisons) {
  de_file <- paste0("results/tables/DE_", comp, ".csv")
  if (file.exists(de_file)) {
    de_res <- read.csv(de_file)
    sig_up <- sum(!is.na(de_res$padj) & de_res$padj < 0.05 & de_res$log2FoldChange > 1, na.rm = TRUE)
    sig_down <- sum(!is.na(de_res$padj) & de_res$padj < 0.05 & de_res$log2FoldChange < -1, na.rm = TRUE)

    de_counts <- rbind(de_counts, data.frame(
      comparison = comp,
      up = sig_up,
      down = sig_down
    ))
  }
}

# Get top 3 upregulated genes from 30m vs 6m (padj < 0.05 & log2FC > 1)
up_genes <- de_results[!is.na(de_results$padj) & de_results$padj < 0.05 & de_results$log2FoldChange > 1, ]
top_up <- head(up_genes[order(up_genes$padj), ], 3)
if (nrow(top_up) == 0) {
  top_up_genes <- "None"
} else {
  top_up_genes <- paste(top_up$gene_symbol, collapse = ", ")
}

# Get top 3 downregulated genes from 30m vs 6m (padj < 0.05 & log2FC < -1)
down_genes <- de_results[!is.na(de_results$padj) & de_results$padj < 0.05 & de_results$log2FoldChange < -1, ]
top_down <- head(down_genes[order(down_genes$padj), ], 3)
if (nrow(top_down) == 0) {
  top_down_genes <- "None"
} else {
  top_down_genes <- paste(top_down$gene_symbol, collapse = ", ")
}

# Get enrichment counts
enrich_file <- "results/tables/gProfiler_GO_KEGG.csv"
if (file.exists(enrich_file)) {
  enrich_data <- read.csv(enrich_file)
  go_terms <- sum(enrich_data$source == "GO:BP", na.rm = TRUE)
  kegg_terms <- sum(enrich_data$source == "KEGG", na.rm = TRUE)
} else {
  go_terms <- 0
  kegg_terms <- 0
}

# Create summary data frame
summary_df <- data.frame(
  Metric = c(
    "Total_samples",
    "Age_groups",
    "Samples_per_group",
    "Total_genes_before_filtering",
    "Total_genes_after_filtering",
    "DE_genes_15m_vs_6m_up",
    "DE_genes_15m_vs_6m_down",
    "DE_genes_24m_vs_15m_up",
    "DE_genes_24m_vs_15m_down",
    "DE_genes_30m_vs_24m_up",
    "DE_genes_30m_vs_24m_down",
    "DE_genes_30m_vs_6m_up",
    "DE_genes_30m_vs_6m_down",
    "Enrichment_GO_terms_up",
    "Enrichment_KEGG_pathways_up",
    "Top_3_upregulated_genes",
    "Top_3_downregulated_genes"
  ),
  Value = c(
    nrow(metadata),
    "6m, 15m, 24m, 30m",
    paste(sprintf("%s=%d", names(samples_per_group), samples_per_group), collapse = "; "),
    format(nrow(counts_raw), big.mark = ","),
    format(nrow(counts_filtered), big.mark = ","),
    de_counts$up[de_counts$comparison == "15m_vs_6m"],
    de_counts$down[de_counts$comparison == "15m_vs_6m"],
    de_counts$up[de_counts$comparison == "24m_vs_15m"],
    de_counts$down[de_counts$comparison == "24m_vs_15m"],
    de_counts$up[de_counts$comparison == "30m_vs_24m"],
    de_counts$down[de_counts$comparison == "30m_vs_24m"],
    de_counts$up[de_counts$comparison == "30m_vs_6m"],
    de_counts$down[de_counts$comparison == "30m_vs_6m"],
    go_terms,
    kegg_terms,
    top_up_genes,
    top_down_genes
  )
)

# Save summary table
write.csv(summary_df, "results/tables/analysis_summary.csv", row.names = FALSE)
message("  Saved: results/tables/analysis_summary.csv")

# Print summary
message("\n  Analysis Summary:")
message(paste(rep("-", 50), collapse = ""))
for (i in 1:nrow(summary_df)) {
  message(sprintf("  %-30s | %s", summary_df$Metric[i], summary_df$Value[i]))
}

# ==============================================================================
# 3. CREATE RESULTS INDEX
# ==============================================================================
message("\nStep 3: Creating results index...")

# Get list of all result files
figures <- list.files("results/figures/", pattern = "\\.png$", full.names = FALSE)
tables <- list.files("results/tables/", pattern = "\\.csv$", full.names = FALSE)

# Create RESULTS_INDEX.md
index_content <- c(
  "# GSE225576 Mouse Brain Aging RNA-seq - Results Index",
  "",
  "## Project Overview",
  "",
  "This analysis examines differential gene expression in mouse brain tissue across four age groups:",
  "- 6 months (young adult)",
  "- 15 months (middle age)",
  "- 24 months (early aged)",
  "- 30 months (aged)",
  "",
  "## Key Findings",
  "",
  paste0("- **Total samples analyzed:** ", nrow(metadata)),
  paste0("- **Genes after filtering:** ", format(nrow(counts_filtered), big.mark = ",")),
  paste0("- **Differential expression (30m vs 6m):** ", de_counts$up[de_counts$comparison == "30m_vs_6m"], " up, ", de_counts$down[de_counts$comparison == "30m_vs_6m"], " down"),
  "- **Major biological themes:**",
  "  - Immune response and inflammation (strongest enrichment)",
  "  - Complement system activation",
  "  - Glial cell activation (astrogliosis)",
  "  - Synaptic transmission downregulation",
  "",
  "## Generated Figures",
  ""
)

# Add figure descriptions
figure_descriptions <- list(
  "qc_library_size.png" = "Library size distribution across samples",
  "qc_pca.png" = "PCA plot showing sample clustering by age group",
  "qc_correlation.png" = "Sample correlation heatmap",
  "volcano_30m_vs_6m.png" = "Volcano plot: 30m vs 6m differential expression",
  "volcano_30m_vs_6m_polished.png" = "Volcano plot with ggrepel labels (polished)",
  "heatmap_top50.png" = "Heatmap of top 50 differentially expressed genes",
  "heatmap_top50_annotated.png" = "Heatmap with gene functional cluster annotations",
  "marker_genes.png" = "Expression trends for aging marker genes",
  "marker_genes_significance.png" = "Marker genes with statistical significance brackets",
  "DE_summary.png" = "Summary of differential expression across all comparisons",
  "enrichment_dotplot.png" = "GO and KEGG enrichment analysis results"
)

for (fig in figures) {
  desc <- figure_descriptions[[fig]]
  if (is.null(desc)) desc <- "Analysis figure"
  index_content <- c(index_content, paste0("- `", fig, "`: ", desc))
}

index_content <- c(index_content, "", "## Generated Tables", "")

# Add table descriptions
table_descriptions <- list(
  "analysis_summary.csv" = "Comprehensive analysis summary metrics",
  "brain_counts.rds" = "Filtered count matrix (processed)",
  "brain_counts_filtered.rds" = "Filtered count matrix (CPM filtered)",
  "brain_metadata.rds" = "Sample metadata with age groups",
  "brain_vst.rds" = "VST-transformed normalized counts",
  "DE_15m_vs_6m.csv" = "Differential expression results: 15m vs 6m",
  "DE_24m_vs_15m.csv" = "Differential expression results: 24m vs 15m",
  "DE_30m_vs_24m.csv" = "Differential expression results: 30m vs 24m",
  "DE_30m_vs_6m.csv" = "Differential expression results: 30m vs 6m",
  "gProfiler_GO_KEGG.csv" = "Functional enrichment results (GO + KEGG)"
)

for (tbl in tables) {
  desc <- table_descriptions[[tbl]]
  if (is.null(desc)) desc <- "Data table"
  index_content <- c(index_content, paste0("- `", tbl, "`: ", desc))
}

index_content <- c(index_content, "", "## Methodology", "",
  "1. **Quality Control:** Library size, PCA, and correlation analysis",
  "2. **Filtering:** Genes with CPM > 1 in at least 3 samples retained",
  "3. **Normalization:** DESeq2 variance stabilizing transformation (VST)",
  "4. **Differential Expression:** DESeq2 with ashr LFC shrinkage",
  "5. **Enrichment:** gProfiler2 (GO:BP + KEGG, FDR corrected)",
  "",
  "## Session Information",
  "",
  "See `session_info.txt` for detailed R environment information."
)

# Write index file
writeLines(index_content, "results/RESULTS_INDEX.md")
message("  Saved: results/RESULTS_INDEX.md")

# ==============================================================================
# 4. UPDATE README.md
# ==============================================================================
message("\nStep 4: Checking for README.md...")

if (file.exists("README.md")) {
  message("  README.md found, appending results summary...")

  readme_appendix <- c(
    "",
    "## Results Summary",
    "",
    "### Overview",
    paste0("Analysis of ", nrow(metadata), " mouse brain samples across 4 age groups (6m, 15m, 24m, 30m) identified ",
           de_counts$up[de_counts$comparison == "30m_vs_6m"] + de_counts$down[de_counts$comparison == "30m_vs_6m"],
           " differentially expressed genes between aged (30m) and young (6m) mice."),
    "",
    "### Key Biological Findings",
    "",
    "- **Neuroinflammation:** Strong upregulation of immune response genes (Clec7a, Gfap, C4b)",
    "- **Complement Activation:** Classical complement pathway components significantly elevated",
    "- **Synaptic Decline:** Downregulation of synaptic transmission markers (Syp, Snap25, Camk2a)",
    "- **Astrocyte Activation:** Gfap and other astrocyte markers show significant age-related increase",
    "",
    "### Main Result Files",
    "",
    "| File | Description |",
    "|------|-------------|",
    "| `results/RESULTS_INDEX.md` | Complete results documentation |",
    "| `results/tables/analysis_summary.csv` | Analysis metrics summary |",
    "| `results/tables/DE_30m_vs_6m.csv` | Full differential expression results |",
    "| `results/figures/volcano_30m_vs_6m_polished.png` | Main volcano plot |",
    "| `results/figures/heatmap_top50_annotated.png` | Gene expression heatmap |",
    "",
    "See `results/RESULTS_INDEX.md` for complete file listing and detailed descriptions."
  )

  # Read existing content and append
  existing_content <- readLines("README.md")
  new_content <- c(existing_content, readme_appendix)
  writeLines(new_content, "README.md")
  message("  Results summary appended to README.md")
} else {
  message("  README.md not found, skipping update")
}

# ==============================================================================
# Complete
# ==============================================================================
message("\n========================================")
message("Session summary generation complete!")
message("========================================")
message("\nGenerated files:")
message("  - results/session_info.txt")
message("  - results/tables/analysis_summary.csv")
message("  - results/RESULTS_INDEX.md")
if (file.exists("README.md")) {
  message("  - README.md (updated)")
}
