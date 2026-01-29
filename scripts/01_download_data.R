#!/usr/bin/env Rscript
# Download and process GSE225576 brain samples
# GSE225576 - Mouse brain aging RNA-seq

message("Loading GEOquery...")
library(GEOquery)

message("Downloading GSE225576 metadata...")
gse <- getGEO("GSE225576", GSEMatrix = TRUE)
eset <- gse[[1]]

message("Extracting metadata...")
metadata <- pData(eset)

# Filter for Brain samples only
message("Filtering for Brain tissue...")
brain_samples <- metadata[metadata$source_name_ch1 == "Brain", ]
message("  Brain samples found: ", nrow(brain_samples))

# Extract age from title (e.g., "Brain_6mo_1" -> "6m")
message("Parsing age groups from titles...")
ages <- gsub("Brain_([0-9]+)mo_.*", "\\1m", brain_samples$title)

# Create count_file_id to match column names in count matrix
# "Brain_6mo_1" -> "6M1"
count_file_id <- gsub("Brain_([0-9]+)mo_([0-9]+)", "\\1M\\2", brain_samples$title)

metadata_clean <- data.frame(
  sample_id = make.names(brain_samples$title, unique = TRUE),
  age_group = factor(ages,
                     levels = c("6m", "15m", "24m", "30m")),
  title = brain_samples$title,
  count_file_id = count_file_id,
  geo_accession = brain_samples$geo_accession
)

message("Age group distribution:")
print(table(metadata_clean$age_group))

# Download supplementary count file
message("\nDownloading supplementary count data...")
supp_dir <- tempdir()
supp_files <- getGEOSuppFiles("GSE225576", baseDir = supp_dir)

# Find and read the Brain count matrix file
count_file <- rownames(supp_files)[grepl("Brain", rownames(supp_files), ignore.case = TRUE)]
if (length(count_file) == 0) {
  stop("Brain count file not found in supplementary files")
}
count_file <- count_file[1]

message("Reading count matrix: ", basename(count_file))
counts <- read.delim(count_file, row.names = 1, check.names = FALSE)
message("  Raw count matrix dimensions: ", nrow(counts), " x ", ncol(counts))

# Filter and reorder count matrix to match metadata
# Count file columns are like "6M1", "6M2" - use count_file_id to match
brain_counts <- counts[, metadata_clean$count_file_id, drop = FALSE]
colnames(brain_counts) <- metadata_clean$sample_id

# Checks
message("\n---")
message("Number of samples: ", ncol(brain_counts))
message("Number of genes: ", nrow(brain_counts))
message("---")
message("Metadata (first rows):")
print(head(metadata_clean))

# Save
message("\nSaving to data/processed/...")
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}
saveRDS(brain_counts, "data/processed/brain_counts.rds")
saveRDS(metadata_clean, "data/processed/brain_metadata.rds")

message("Done.")
