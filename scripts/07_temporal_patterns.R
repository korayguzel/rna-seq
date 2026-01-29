#!/usr/bin/env Rscript
# Temporal pattern analysis using maSigPro
# GSE225576 - Mouse brain aging RNA-seq

message("Loading packages...")
library(tidyverse)

# Check/install maSigPro
if (!requireNamespace("maSigPro", quietly = TRUE)) {
  message("Installing maSigPro...")
  BiocManager::install("maSigPro", ask = FALSE, update = FALSE)
}
library(maSigPro)
library(pheatmap)

# Set seed for reproducibility
set.seed(123)

# Create output directories
if (!dir.exists("results/tables")) {
  dir.create("results/tables", recursive = TRUE)
}
if (!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================
message("\nStep 1: Loading and preparing data...")

# Load count matrix and metadata
counts <- readRDS("data/processed/brain_counts_filtered.rds")
metadata <- readRDS("data/processed/brain_metadata.rds")

message("  Counts: ", nrow(counts), " genes x ", ncol(counts), " samples")

# Ensure counts are integers (required by maSigPro)
if (!all(counts == floor(counts), na.rm = TRUE)) {
  message("  Rounding counts to integers...")
  counts <- round(counts)
}

# Create Edesign (experimental design matrix)
# Time: numeric age in months
# Replicate: sample identifier
# Group: all samples belong to same group (single timecourse)

# Map age groups to numeric time
age_map <- c("6m" = 6, "15m" = 15, "24m" = 24, "30m" = 30)
metadata$Time <- age_map[as.character(metadata$age_group)]

# Create replicate IDs based on age group (4 replicates per timepoint)
metadata$Replicate <- ave(1:nrow(metadata), metadata$age_group, FUN = seq_along)

# For maSigPro single timecourse, we use Group as time indicator
# Set all to 1 for single timecourse analysis
metadata$Group <- rep(1, nrow(metadata))

# Create Edesign matrix
edesign <- data.frame(
  Time = metadata$Time,
  Replicate = metadata$Replicate,
  Group = metadata$Group,
  row.names = metadata$sample_id
)

# Reorder counts to match edesign row order
counts <- counts[, rownames(edesign), drop = FALSE]

message("  Experimental design:")
print(table(edesign$Time))

# Create design matrix for maSigPro (polynomial degree 2 - appropriate for 4 timepoints)
message("\n  Creating design matrix (degree = 2)...")
design <- make.design.matrix(edesign, degree = 2)
message("  Design matrix created with ", ncol(design$groups), " columns")
message("  Design matrix dimensions: ", nrow(design$groups), " x ", ncol(design$groups))

# ==============================================================================
# 2. MODEL FITTING
# ==============================================================================
message("\nStep 2: Fitting temporal models...")

# Step 1: Identify genes with significant temporal changes
message("  Running p.vector (Q = 0.05)...")
pvector_result <- p.vector(
  counts,
  design,
  Q = 0.05,
  MT.adjust = "BH",
  min.obs = 3
)

sig_genes_count <- nrow(pvector_result$SELEC)
message("  Genes with significant temporal patterns: ", sig_genes_count)

if (sig_genes_count < 10) {
  message("\n  WARNING: Very few genes passed significance threshold.")
  message("  Consider increasing Q threshold or interpreting results as exploratory.")
}

# Step 2: Fit polynomial models to significant genes
# Note: T.fit may fail with single-group design; handle gracefully
message("\n  Running T.fit (rsq = 0.6)...")
tfit_result <- tryCatch({
  T.fit(
    pvector_result,
    design,
    step.method = "backward",
    min.obs = 3
  )
}, error = function(e) {
  message("  T.fit failed: ", conditionMessage(e))
  message("  Using p.vector results directly for clustering...")
  return(NULL)
})

# Step 3: Extract significant genes
if (!is.null(tfit_result)) {
  message("\n  Extracting significant genes (rsq = 0.6)...")
  sig_genes <- get.siggenes(
    tfit_result,
    rsq = 0.6,
    vars = "groups"
  )
  final_sig_count <- length(sig_genes$sig.genes$`1`$g)
} else {
  # Fallback: use p.vector selected genes directly
  sig_genes <- NULL
  final_sig_count <- sig_genes_count
}
message("  Final significant genes: ", final_sig_count)

# ==============================================================================
# 3. CLUSTERING
# ==============================================================================
message("\nStep 3: Clustering temporal patterns...")

# Extract fitted values for clustering
if (final_sig_count > 0) {

  # Get gene names from appropriate source
  if (!is.null(sig_genes)) {
    sig_gene_names <- sig_genes$sig.genes$`1`$g
  } else {
    # Use p.vector selected genes
    sig_gene_names <- rownames(pvector_result$SELEC)
  }

  # Get fitted values or use expression profiles
  if (!is.null(sig_genes) && !is.null(sig_genes$sig.genes$`1`$sig.profiles)) {
    profiles <- sig_genes$sig.genes$`1`$sig.profiles
  } else if (!is.null(tfit_result) && !is.null(tfit_result$coefficients)) {
    # Extract from tfit_result if available
    profiles <- tfit_result$coefficients[sig_gene_names, , drop = FALSE]
  } else {
    # Use original expression values
    profiles <- NULL
  }

  # Determine number of clusters
  k <- min(5, final_sig_count)
  if (k < 2) {
    k <- 2
  }

  message("  Clustering into ", k, " groups...")

  # Perform clustering on fitted values (or use scaled profiles)
  if (is.null(profiles) || nrow(profiles) == 0) {
    # Fallback: use original expression values
    profiles <- counts[sig_gene_names, , drop = FALSE]
    # Average by timepoint
    profile_means <- data.frame(
      gene = sig_gene_names,
      Time6 = rowMeans(profiles[, metadata$Time == 6, drop = FALSE]),
      Time15 = rowMeans(profiles[, metadata$Time == 15, drop = FALSE]),
      Time24 = rowMeans(profiles[, metadata$Time == 24, drop = FALSE]),
      Time30 = rowMeans(profiles[, metadata$Time == 30, drop = FALSE])
    )
    rownames(profile_means) <- sig_gene_names
    profiles <- as.matrix(profile_means[, -1])
  }

  # Scale profiles for clustering
  profiles_scaled <- t(scale(t(profiles)))
  profiles_scaled <- profiles_scaled[complete.cases(profiles_scaled), ]

  # K-means clustering
  kmeans_result <- kmeans(profiles_scaled, centers = k, nstart = 25)
  clusters <- kmeans_result$cluster

  message("  Cluster assignments:")
  for (i in 1:k) {
    cluster_genes <- names(clusters)[clusters == i]
    message("    Cluster ", i, ": ", length(cluster_genes), " genes")
  }

  # ==============================================================================
  # 4. BIOLOGICAL INTERPRETATION
  # ==============================================================================
  message("\nStep 4: Biological interpretation...")

  # Identify cluster profiles
  cluster_profiles <- data.frame()
  for (i in 1:k) {
    cluster_genes <- names(clusters)[clusters == i]
    if (length(cluster_genes) == 0) next

    # Get mean profile for cluster
    cluster_data <- profiles_scaled[cluster_genes, , drop = FALSE]
    mean_profile <- colMeans(cluster_data, na.rm = TRUE)

    # Determine pattern type
    pattern <- "Complex"
    if (all(diff(mean_profile) > -0.1)) {
      pattern <- "Linear Increase"
    } else if (all(diff(mean_profile) < 0.1)) {
      pattern <- "Linear Decrease"
    } else if (which.max(mean_profile) <= 2) {
      pattern <- "Early Response"
    } else if (which.max(mean_profile) >= 4) {
      pattern <- "Late Response"
    } else if (which.max(mean_profile) %in% c(2, 3)) {
      pattern <- "Transient"
    }

    cluster_profiles <- rbind(cluster_profiles, data.frame(
      cluster = i,
      pattern = pattern,
      n_genes = length(cluster_genes),
      Time6 = mean_profile[1],
      Time15 = mean_profile[2],
      Time24 = mean_profile[3],
      Time30 = mean_profile[4]
    ))

    message("\n  Cluster ", i, " (", pattern, "): ", length(cluster_genes), " genes")

    # Top 5 genes (by cluster assignment order if rsq not available)
    if (!is.null(tfit_result) && !is.null(tfit_result$rsq)) {
      rsq_values <- tfit_result$rsq[sig_gene_names]
      cluster_rsq <- rsq_values[cluster_genes]
      top5 <- head(sort(cluster_rsq, decreasing = TRUE), 5)
      message("    Top 5 genes: ", paste(names(top5), collapse = ", "))
    } else {
      message("    Genes: ", paste(head(cluster_genes, 5), collapse = ", "),
              ifelse(length(cluster_genes) > 5, "...", ""))
    }
  }

  # Check for expected markers
  message("\n  Checking expected marker patterns...")
  up_markers <- c("Gfap", "C4b", "Ctss", "Clec7a")
  down_markers <- c("Syp", "Snap25", "Camk2a")

  for (marker in up_markers) {
    if (marker %in% sig_gene_names) {
      cl <- clusters[marker]
      pattern <- cluster_profiles$pattern[cluster_profiles$cluster == cl]
      message("    ", marker, ": Cluster ", cl, " (", pattern, ")")
    } else if (marker %in% rownames(counts)) {
      message("    ", marker, ": Not significant")
    } else {
      message("    ", marker, ": Not in dataset")
    }
  }

  for (marker in down_markers) {
    if (marker %in% sig_gene_names) {
      cl <- clusters[marker]
      pattern <- cluster_profiles$pattern[cluster_profiles$cluster == cl]
      message("    ", marker, ": Cluster ", cl, " (", pattern, ")")
    } else if (marker %in% rownames(counts)) {
      message("    ", marker, ": Not significant")
    } else {
      message("    ", marker, ": Not in dataset")
    }
  }

  # ==============================================================================
  # 5. SAVE OUTPUTS
  # ==============================================================================
  message("\nStep 5: Saving outputs...")

  # Save significant genes
  if (!is.null(tfit_result) && !is.null(tfit_result$rsq)) {
    rsq_vals <- tfit_result$rsq[sig_gene_names]
  } else {
    rsq_vals <- rep(NA, length(sig_gene_names))
  }
  sig_genes_df <- data.frame(
    gene_id = sig_gene_names,
    rsq = rsq_vals,
    cluster = clusters[sig_gene_names]
  )
  write.csv(sig_genes_df, "results/tables/maSigPro_significant_genes.csv", row.names = FALSE)
  message("  Saved: results/tables/maSigPro_significant_genes.csv")

  # Save cluster mapping
  cluster_df <- data.frame(
    gene_id = names(clusters),
    cluster = clusters
  )
  write.csv(cluster_df, "results/tables/maSigPro_clusters.csv", row.names = FALSE)
  message("  Saved: results/tables/maSigPro_clusters.csv")

  # Save cluster profiles
  write.csv(cluster_profiles, "results/tables/cluster_profiles.csv", row.names = FALSE)
  message("  Saved: results/tables/cluster_profiles.csv")

  # Create cluster visualization
  message("\n  Generating cluster plot...")

  # Prepare data for plotting
  plot_data <- data.frame()
  for (i in 1:k) {
    cluster_genes <- names(clusters)[clusters == i]
    if (length(cluster_genes) == 0) next

    cluster_data <- profiles_scaled[cluster_genes, , drop = FALSE]
    mean_profile <- colMeans(cluster_data, na.rm = TRUE)
    se_profile <- apply(cluster_data, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))

    cluster_plot_data <- data.frame(
      Time = c(6, 15, 24, 30),
      Mean = mean_profile,
      SE = se_profile,
      Cluster = paste0("Cluster ", i, " (", cluster_profiles$pattern[i], ")"),
      ClusterNum = i
    )
    plot_data <- rbind(plot_data, cluster_plot_data)
  }

  p_clusters <- ggplot(plot_data, aes(x = Time, y = Mean, color = Cluster, fill = Cluster)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE), alpha = 0.2, color = NA) +
    scale_x_continuous(breaks = c(6, 15, 24, 30)) +
    labs(x = "Age (months)", y = "Standardized Expression",
         title = "Temporal Expression Patterns",
         subtitle = paste0("maSigPro clustering (k=", k, ")")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")

  ggsave("results/figures/maSigPro_clusters.png", p_clusters,
         dpi = 300, width = 10, height = 6, units = "in", bg = "white")
  message("  Saved: results/figures/maSigPro_clusters.png")

  # Create heatmap of significant genes
  message("\n  Generating heatmap...")

  # Order genes by cluster
  gene_order <- names(clusters)[order(clusters)]
  heatmap_data <- profiles_scaled[gene_order, , drop = FALSE]

  # Create cluster annotation
  annotation_row <- data.frame(
    Cluster = factor(clusters[gene_order]),
    row.names = gene_order
  )

  # Create timepoint annotation
  annotation_col <- data.frame(
    Age = factor(c("6m", "15m", "24m", "30m")),
    row.names = colnames(heatmap_data)
  )

  png("results/figures/maSigPro_heatmap.png", width = 8, height = 10, units = "in", res = 300)
  pheatmap(heatmap_data,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_row = annotation_row,
           annotation_col = annotation_col,
           show_rownames = FALSE,
           fontsize_col = 12,
           main = "Temporal Patterns - Significant Genes",
           color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
  message("  Saved: results/figures/maSigPro_heatmap.png")

} else {
  message("\n  WARNING: No genes passed rsq threshold.")
  message("  Try lowering rsq threshold or check data quality.")
}

# ==============================================================================
# Complete
# ==============================================================================
message("\n========================================")
message("Temporal pattern analysis complete!")
message("========================================")
message("\nNote: Small sample size (n=16) makes results exploratory.")
message("Consider validating patterns with additional replicates.")
