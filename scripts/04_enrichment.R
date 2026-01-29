#!/usr/bin/env Rscript
# Functional enrichment analysis using gProfiler2
# GSE225576 - Mouse brain aging RNA-seq (30m vs 6m comparison)

message("Loading packages...")
library(gprofiler2)
library(tidyverse)
library(patchwork)

# Create directories if needed
if (!dir.exists("results/tables")) {
  dir.create("results/tables", recursive = TRUE)
}
if (!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}

message("Loading DE results...")
de_results <- read.csv("results/tables/DE_30m_vs_6m.csv")

# 1. Prepare upregulated gene list
message("\nPreparing upregulated gene list...")
up_genes <- de_results[!is.na(de_results$padj) &
                         de_results$padj < 0.05 &
                         de_results$log2FoldChange > 1, ]$gene_symbol

message("  Upregulated genes: ", length(up_genes))

if (length(up_genes) < 5) {
  stop("ERROR: Too few upregulated genes for enrichment (< 5)")
}

# 2. Run gProfiler2 enrichment
message("\nRunning gProfiler2 enrichment...")
message("  Organism: mmusculus (mouse)")
message("  Sources: GO:BP, KEGG")
message("  Correction: FDR")

enrich_result <- gost(
  query = up_genes,
  organism = "mmusculus",
  sources = c("GO:BP", "KEGG"),
  correction_method = "fdr",
  significant = TRUE,
  user_threshold = 0.05
)

# Check if results exist
if (is.null(enrich_result) || is.null(enrich_result$result) || nrow(enrich_result$result) == 0) {
  message("\n  WARNING: No enriched terms found!")
  message("  This is unexpected - upregulated genes may not map to known pathways.")
  quit(status = 0)
}

results_df <- enrich_result$result
message("  Enriched terms found: ", nrow(results_df))

# 3. Print summary
message("\n", paste(rep("=", 60), collapse = ""))
message("ENRICHMENT SUMMARY")
message(paste(rep("=", 60), collapse = ""))

# Count by source
go_terms <- results_df[results_df$source == "GO:BP", ]
kegg_terms <- results_df[results_df$source == "KEGG", ]

message("GO Biological Process terms: ", nrow(go_terms))
message("KEGG pathways: ", nrow(kegg_terms))

# Top enriched terms (by lowest p-value)
message("\nTop 10 enriched terms:")
top_terms <- head(results_df[order(results_df$p_value), ], 10)
for (i in 1:nrow(top_terms)) {
  term <- top_terms[i, ]
  source_label <- ifelse(term$source == "GO:BP", "[GO]", "[KEGG]")
  message(sprintf("  %d. %s %s (p=%.2e, n=%d)",
                  i, source_label, term$term_name,
                  term$p_value, term$intersection_size))
}

# Check for expected terms (neuroinflammation, complement, etc.)
message("\nChecking for expected aging terms...")
expected_patterns <- c("inflammatory", "immune", "complement", "glia",
                       "neuroinflammation", "cytokine", "macrophage", "microglia")
found_terms <- c()
for (pattern in expected_patterns) {
  matches <- results_df[grep(pattern, results_df$term_name, ignore.case = TRUE), ]
  if (nrow(matches) > 0) {
    found_terms <- c(found_terms, paste0(pattern, " (", nrow(matches), ")"))
  }
}
if (length(found_terms) > 0) {
  message("  Found terms related to: ", paste(found_terms, collapse = ", "))
} else {
  message("  Note: No direct matches for expected aging terms found")
}

# 4. Save results
output_csv <- "results/tables/gProfiler_GO_KEGG.csv"

# Convert list columns to character for CSV export
results_df_export <- results_df
list_cols <- sapply(results_df_export, is.list)
if (any(list_cols)) {
  for (col in names(list_cols)[list_cols]) {
    results_df_export[[col]] <- sapply(results_df_export[[col]], paste, collapse = ";")
  }
}
write.csv(results_df_export, output_csv, row.names = FALSE)
message("\n  Saved: ", output_csv)

# 5. Create dotplot (top 15 terms)
message("\nCreating dotplot...")

# Select top 15 terms by p-value
top15 <- head(results_df[order(results_df$p_value), ], 15)

# Add source label
top15$source_label <- ifelse(top15$source == "GO:BP", "GO:BP", "KEGG")

# Create plot
p <- ggplot(top15, aes(x = intersection_size,
                       y = reorder(term_name, -p_value),
                       color = p_value,
                       size = intersection_size)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", trans = "log10",
                       name = "p-value") +
  scale_size_continuous(name = "Gene count") +
  facet_grid(source_label ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Intersection Size (Genes)",
       y = "Enriched Term") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 14)
  )

# Add title using patchwork (cleaner than cowplot)
p_with_title <- p + plot_annotation(
  title = "Functional Enrichment - Upregulated Genes (30m vs 6m)",
  theme = theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  )
)

ggsave("results/figures/enrichment_dotplot.png", p_with_title,
       dpi = 300, width = 10, height = 8, units = "in", bg = "white")
message("  Saved: results/figures/enrichment_dotplot.png")

message("\n", paste(rep("=", 60), collapse = ""))
message("Enrichment analysis complete.")
message(paste(rep("=", 60), collapse = ""))
