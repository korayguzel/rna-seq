#!/usr/bin/env Rscript
# Install required packages for mouse brain aging RNA-seq analysis
# GSE225576 - Brain tissue only

message("Starting package installation...")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  message("installed: BiocManager")
}

# CRAN packages
cran_packages <- c("tidyverse", "pheatmap", "ggrepel", "WGCNA")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org/")
      message("installed: ", pkg)
    }, error = function(e) {
      message("FAILED: ", pkg, " - ", conditionMessage(e))
    })
  } else {
    message("installed: ", pkg)
  }
}

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "GEOquery",
  "clusterProfiler",
  "org.Mm.eg.db",
  "EnhancedVolcano",
  "maSigPro"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      message("installed: ", pkg)
    }, error = function(e) {
      message("FAILED: ", pkg, " - ", conditionMessage(e))
    })
  } else {
    message("installed: ", pkg)
  }
}

message("Package installation complete.")
