#!/usr/bin/env Rscript
# Create folder structure for mouse brain aging RNA-seq project

message("Creating folder structure...")

dirs <- c(
  "data/raw",
  "data/processed",
  "scripts",
  "results/figures",
  "results/tables"
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  message("  Created: ", d)
}

message("Done.")
