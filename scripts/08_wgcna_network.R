#!/usr/bin/env Rscript
# WGCNA Co-expression Network Analysis
# GSE225576 - Mouse brain aging RNA-seq
# Note: Small sample size (n=16) makes results exploratory

message("Loading packages...")
library(tidyverse)

# Check/install WGCNA
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  message("Installing WGCNA...")
  install.packages("WGCNA", repos = "https://cloud.r-project.org/")
}
library(WGCNA)
library(pheatmap)

# WGCNA settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

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
message("\nStep 1: Data preparation...")

# Load VST data
vst_matrix <- readRDS("data/processed/brain_vst.rds")
metadata <- readRDS("data/processed/brain_metadata.rds")

message("  Loaded VST data: ", nrow(vst_matrix), " genes x ", ncol(vst_matrix), " samples")

# Transpose: samples as rows, genes as columns (WGCNA standard)
datExpr <- t(vst_matrix)

message("  Transposed: ", nrow(datExpr), " samples x ", ncol(datExpr), " genes")

# Filter top 5000 most variable genes
gene_vars <- apply(datExpr, 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
datExpr <- datExpr[, top_genes]

message("  Selected top ", ncol(datExpr), " variable genes")

# Check for missing values
if (any(is.na(datExpr))) {
  message("  WARNING: Missing values detected. Removing genes with NAs...")
  datExpr <- datExpr[, colSums(is.na(datExpr)) == 0]
  message("  Retained ", ncol(datExpr), " genes after NA removal")
}

# Check for outlier samples
message("\n  Checking for outlier samples...")
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot sample dendrogram
png("results/figures/wgcna_dendrogram.png", width = 10, height = 6, units = "in", res = 300)
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     xlab = "", sub = "", cex.lab = 1.5, cex.axis = 1.5)
abline(h = 15, col = "red", lty = 2)
dev.off()
message("  Saved: results/figures/wgcna_dendrogram.png")

# Identify outliers (height > 15)
outlier_cutoff <- 15
outliers <- sampleTree$labels[sampleTree$height > outlier_cutoff]
if (length(outliers) > 0) {
  message("  Detected ", length(outliers), " outlier samples: ", paste(outliers, collapse = ", "))
  message("  Removing outliers...")
  datExpr <- datExpr[!(rownames(datExpr) %in% outliers), ]
} else {
  message("  No outliers detected")
}

message("\n  Final data for network construction:")
message("    Samples: ", nrow(datExpr))
message("    Genes: ", ncol(datExpr))

# ==============================================================================
# 2. SOFT THRESHOLD SELECTION
# ==============================================================================
message("\nStep 2: Selecting soft thresholding power...")

# Choose powers to test
powers <- c(1:20)

# Call network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot results
png("results/figures/wgcna_soft_threshold.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))

# Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale-free topology fit", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.8, col = "red", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

dev.off()
message("  Saved: results/figures/wgcna_soft_threshold.png")

# Select power where R^2 > 0.8 (or closest available)
r2_values <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
valid_powers <- sft$fitIndices[r2_values > 0.8, 1]
if (length(valid_powers) > 0) {
  softPower <- min(valid_powers)
} else {
  softPower <- 6  # Default if no power reaches 0.8
  message("  WARNING: No power reached R^2 > 0.8. Using default power = 6")
}

message("  Selected soft threshold power: ", softPower)
message("  Scale-free topology R^2 at selected power: ",
        round(r2_values[sft$fitIndices[, 1] == softPower], 3))

# ==============================================================================
# 3. NETWORK CONSTRUCTION
# ==============================================================================
message("\nStep 3: Constructing network...")
message("  This may take 5-10 minutes...")

# Build network using blockwiseModules
net <- tryCatch({
  blockwiseModules(
    datExpr,
    power = softPower,
    TOMType = "signed",
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    maxBlockSize = 5000,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    verbose = 3
  )
}, error = function(e) {
  message("  Error in blockwiseModules: ", conditionMessage(e))
  message("  Trying with adjusted parameters...")
  blockwiseModules(
    datExpr,
    power = softPower,
    TOMType = "unsigned",
    minModuleSize = 20,
    mergeCutHeight = 0.30,
    maxBlockSize = 5000,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    verbose = 3
  )
})

# Module information
module_labels <- net$colors
module_colors <- labels2colors(module_labels)
n_modules <- length(unique(module_labels)) - 1  # Exclude grey (unassigned)

message("\n  Network construction complete!")
message("  Number of modules identified: ", n_modules)

# Module sizes
table_colors <- table(module_colors)
message("  Module sizes:")
for (color in names(table_colors)) {
  if (color != "grey") {
    message("    ", color, ": ", table_colors[color], " genes")
  }
}

# ==============================================================================
# 4. MODULE-AGE CORRELATION
# ==============================================================================
message("\nStep 4: Correlating modules with age...")

# Create trait matrix (Age as continuous variable)
age_map <- c("6m" = 6, "15m" = 15, "24m" = 24, "30m" = 30)
sample_ages <- age_map[as.character(metadata$age_group[match(rownames(datExpr), metadata$sample_id)])]

# Create trait dataframe
traitData <- data.frame(Age = sample_ages)
rownames(traitData) <- rownames(datExpr)

# Calculate module eigengenes
MEs <- net$MEs

# Rename MEs to module colors
MEs <- orderMEs(MEs)

# Calculate correlations
moduleTraitCor <- cor(MEs, traitData$Age, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Create correlation dataframe
cor_df <- data.frame(
  Module = rownames(moduleTraitCor),
  Age_Correlation = as.vector(moduleTraitCor),
  P_value = as.vector(moduleTraitPvalue)
)
cor_df$Significant <- abs(cor_df$Age_Correlation) > 0.5 & cor_df$P_value < 0.05

write.csv(cor_df, "results/tables/wgcna_module_trait_cor.csv", row.names = FALSE)
message("  Saved: results/tables/wgcna_module_trait_cor.csv")

# Print significant correlations
message("\n  Module-age correlations:")
sig_modules <- cor_df[cor_df$Significant, ]
if (nrow(sig_modules) > 0) {
  for (i in 1:nrow(sig_modules)) {
    direction <- ifelse(sig_modules$Age_Correlation[i] > 0, "positive", "negative")
    message("    ", sig_modules$Module[i], ": ",
            round(sig_modules$Age_Correlation[i], 3),
            " (", direction, ", p = ",
            format(sig_modules$P_value[i], digits = 2), ")")
  }
} else {
  message("    No modules significantly correlated with age at |cor| > 0.5, p < 0.05")
  # Print top correlations anyway
  top_cor <- cor_df[order(abs(cor_df$Age_Correlation), decreasing = TRUE), ][1:3, ]
  message("    Top 3 correlations:")
  for (i in 1:nrow(top_cor)) {
    message("      ", top_cor$Module[i], ": ", round(top_cor$Age_Correlation[i], 3))
  }
}

# Plot module-trait correlation heatmap
png("results/figures/wgcna_module_trait.png", width = 6, height = 8, units = "in", res = 300)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8, 4, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Age",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-Age Correlations")

dev.off()
message("  Saved: results/figures/wgcna_module_trait.png")

# Plot eigengene heatmap
png("results/figures/wgcna_eigengene_heatmap.png", width = 10, height = 8, units = "in", res = 300)

# Order samples by age
sample_order <- order(traitData$Age)
MEs_ordered <- MEs[sample_order, ]

# Create age annotation
age_annot <- data.frame(
  Age = factor(metadata$age_group[match(rownames(MEs_ordered), metadata$sample_id)]),
  row.names = rownames(MEs_ordered)
)

# Convert MEs to matrix for pheatmap
ME_matrix <- as.matrix(MEs_ordered)

pheatmap(t(ME_matrix),
         annotation_col = age_annot,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         main = "Module Eigengenes by Sample",
         color = colorRampPalette(c("blue", "white", "red"))(100))

dev.off()
message("  Saved: results/figures/wgcna_eigengene_heatmap.png")

# ==============================================================================
# 5. HUB GENE IDENTIFICATION
# ==============================================================================
message("\nStep 5: Identifying hub genes...")

# Get gene-module assignments
gene_names <- colnames(datExpr)
gene_module_df <- data.frame(
  gene_id = gene_names,
  module_color = module_colors,
  stringsAsFactors = FALSE
)

# Calculate intramodular connectivity
message("  Calculating intramodular connectivity...")
connectivity <- intramodularConnectivity.fromExpr(datExpr, module_colors)
gene_module_df$kWithin <- connectivity$kWithin[match(gene_names, rownames(connectivity))]

# Calculate gene significance (correlation with age)
message("  Calculating gene significance...")
geneTraitSignificance <- as.data.frame(cor(datExpr, traitData$Age, use = "p"))
gene_module_df$GS <- abs(geneTraitSignificance[match(gene_names, rownames(geneTraitSignificance)), 1])

# Calculate module membership
message("  Calculating module membership...")
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))

# Match gene to its module's membership
get_MM <- function(gene, module) {
  if (module == "grey") return(NA)
  mod_col <- paste0("ME", module)
  if (mod_col %in% names(geneModuleMembership)) {
    return(abs(geneModuleMembership[gene, mod_col]))
  }
  return(NA)
}

gene_module_df$MM <- mapply(get_MM, gene_module_df$gene_id, gene_module_df$module_color)

# Identify hub genes for significant modules
hub_genes_df <- data.frame()

if (nrow(sig_modules) > 0) {
  for (mod in sig_modules$Module) {
    mod_color <- substring(mod, 3)
    mod_genes <- gene_module_df[gene_module_df$module_color == mod_color, ]

    if (nrow(mod_genes) == 0) next

    # Hub criteria: MM > 0.8 AND GS > 0.5, OR top 10% connectivity
    conn_cutoff <- quantile(mod_genes$kWithin, 0.9, na.rm = TRUE)

    hubs <- mod_genes[
      (mod_genes$MM > 0.8 & mod_genes$GS > 0.5) |
        mod_genes$kWithin >= conn_cutoff,
    ]

    if (nrow(hubs) > 0) {
      hubs <- hubs[order(hubs$kWithin, decreasing = TRUE), ]
      hubs$rank <- 1:nrow(hubs)
      hubs$hub_criteria <- ifelse(hubs$MM > 0.8 & hubs$GS > 0.5, "MM+GS", "Connectivity")
      hub_genes_df <- rbind(hub_genes_df, hubs)

      message("\n  ", mod, " (", mod_color, ") hub genes:")
      top10 <- head(hubs, 10)
      for (i in 1:nrow(top10)) {
        message("    ", i, ". ", top10$gene_id[i],
                " (MM=", round(top10$MM[i], 2),
                ", GS=", round(top10$GS[i], 2), ")")
      }
    }
  }
} else {
  message("  No significant modules for hub gene analysis")

  # Report top genes in largest module
  largest_module <- names(which.max(table_colors[names(table_colors) != "grey"]))
  mod_genes <- gene_module_df[gene_module_df$module_color == largest_module, ]
  mod_genes <- mod_genes[order(mod_genes$kWithin, decreasing = TRUE), ]

  message("\n  Top 10 genes in largest module (", largest_module, "):")
  for (i in 1:min(10, nrow(mod_genes))) {
    message("    ", i, ". ", mod_genes$gene_id[i],
            " (kWithin=", round(mod_genes$kWithin[i], 2), ")")
  }
}

# Save gene-module assignments
write.csv(gene_module_df, "results/tables/wgcna_module_assignments.csv", row.names = FALSE)
message("\n  Saved: results/tables/wgcna_module_assignments.csv")

# Save hub genes
if (nrow(hub_genes_df) > 0) {
  write.csv(hub_genes_df, "results/tables/wgcna_hub_genes.csv", row.names = FALSE)
  message("  Saved: results/tables/wgcna_hub_genes.csv")
} else {
  message("  No hub genes identified with stringent criteria")
  # Save empty file with headers
  write.csv(data.frame(
    gene_id = character(),
    module_color = character(),
    kWithin = numeric(),
    GS = numeric(),
    MM = numeric(),
    rank = integer(),
    hub_criteria = character()
  ), "results/tables/wgcna_hub_genes.csv", row.names = FALSE)
  message("  Saved empty hub genes file")
}

# ==============================================================================
# 6. BIOLOGICAL INTERPRETATION
# ==============================================================================
message("\nStep 6: Biological interpretation...")

message("\n", paste(rep("=", 60), collapse = ""))
message("WGCNA SUMMARY")
message(paste(rep("=", 60), collapse = ""))

message("\nNetwork Statistics:")
message("  - Samples: ", nrow(datExpr))
message("  - Genes: ", ncol(datExpr))
message("  - Soft threshold power: ", softPower)
message("  - Modules: ", n_modules)

message("\nModule Summary:")
for (color in names(table_colors)) {
  if (color != "grey") {
    cor_val <- cor_df$Age_Correlation[cor_df$Module == paste0("ME", color)]
    if (length(cor_val) > 0) {
      message("  ", color, " (n=", table_colors[color], "): Age cor = ", round(cor_val, 3))
    }
  }
}

# Check for expected patterns
message("\nChecking expected marker genes...")
up_markers <- c("Gfap", "C4b", "Ctss", "Clec7a")
down_markers <- c("Syp", "Snap25", "Camk2a")

for (marker in c(up_markers, down_markers)) {
  if (marker %in% gene_module_df$gene_id) {
    mod <- gene_module_df$module_color[gene_module_df$gene_id == marker]
    cor_val <- cor_df$Age_Correlation[cor_df$Module == paste0("ME", mod)]
    message("  ", marker, ": ", mod, " module (Age cor = ", round(cor_val, 3), ")")
  } else {
    message("  ", marker, ": Not in top 5000 variable genes")
  }
}

message("\n", paste(rep("=", 60), collapse = ""))
message("WARNING: Results are exploratory due to small sample size (n=16)")
message(paste(rep("=", 60), collapse = ""))

# ==============================================================================
# Complete
# ==============================================================================
message("\n========================================")
message("WGCNA analysis complete!")
message("========================================")
message("\nGenerated files:")
message("  Tables:")
message("    - results/tables/wgcna_module_assignments.csv")
message("    - results/tables/wgcna_module_trait_cor.csv")
message("    - results/tables/wgcna_hub_genes.csv")
message("  Figures:")
message("    - results/figures/wgcna_soft_threshold.png")
message("    - results/figures/wgcna_dendrogram.png")
message("    - results/figures/wgcna_module_trait.png")
message("    - results/figures/wgcna_eigengene_heatmap.png")
