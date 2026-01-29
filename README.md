# Mouse Brain Aging Transcriptome Analysis (GSE225576)

A reproducible RNA-Seq analysis pipeline examining transcriptional changes in mouse brain tissue across the aging process. This project demonstrates proficiency in bioinformatics workflows including quality control, differential expression analysis, functional enrichment, and temporal pattern identification.

## Project Overview

This analysis investigates age-related transcriptional changes in mouse brain tissue using bulk RNA-Seq data from the Gene Expression Omnibus (GEO accession: [GSE225576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225576)). The dataset comprises 16 samples across four age groups (6, 15, 24, and 30 months; n=4 per group), enabling characterization of both acute and cumulative aging effects on the brain transcriptome.

**Organism:** *Mus musculus* (mouse)
**Tissue:** Brain
**Platform:** Illumina NovaSeq 6000 (paired-end 2x151 bp)

## Analysis Workflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         RNA-Seq Analysis Pipeline                           │
└─────────────────────────────────────────────────────────────────────────────┘

    ┌──────────────┐
    │  GEO Data    │ ─── GSE225576 (GEOquery)
    │  Retrieval   │
    └──────┬───────┘
           │
           ▼
    ┌──────────────┐     • Library size distribution
    │  Quality     │ ─── • PCA (top 500 variable genes)
    │  Control     │     • Sample correlation heatmap
    └──────┬───────┘
           │
           ▼
    ┌──────────────┐     • CPM > 1 in ≥3 samples
    │  Filtering & │ ─── • 55,285 → 15,893 genes retained
    │  Normalization│    • DESeq2 VST transformation
    └──────┬───────┘
           │
           ├─────────────────────┬─────────────────────┐
           ▼                     ▼                     ▼
    ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
    │  Differential│     │  Temporal    │     │   WGCNA      │
    │  Expression  │     │  Patterns    │     │   Network    │
    │  (DESeq2)    │     │  (maSigPro)  │     │  (exploratory)│
    └──────┬───────┘     └──────┬───────┘     └──────────────┘
           │                     │                    *
           ▼                     ▼             *see Limitations
    ┌──────────────┐     ┌──────────────┐
    │  Functional  │     │   Cluster    │
    │  Enrichment  │     │  Assignment  │
    │  (gProfiler2)│     │  (k-means)   │
    └──────┬───────┘     └──────────────┘
           │
           ▼
    ┌──────────────────────────────────────────────────┐
    │              Publication-Quality Figures          │
    │  • Volcano plots     • Heatmaps                  │
    │  • Enrichment plots  • Marker gene boxplots      │
    └──────────────────────────────────────────────────┘
```

## Workflow and Tools

| Step | Tool/Package | Purpose | Reference |
|------|--------------|---------|-----------|
| Data Retrieval | GEOquery | Download count matrix from GEO | Davis & Meltzer, 2007 |
| Quality Control | ggplot2 | Library size, PCA, correlation plots | Wickham, 2016 |
| Normalization | DESeq2 | Variance stabilizing transformation (VST) | Love *et al.*, 2014 |
| Differential Expression | DESeq2 + ashr | Negative binomial GLM with LFC shrinkage | Love *et al.*, 2014; Stephens, 2017 |
| Functional Enrichment | gProfiler2 | GO:BP and KEGG pathway analysis | Kolberg *et al.*, 2020 |
| Temporal Analysis | maSigPro | Polynomial regression for time-course data | Conesa *et al.*, 2006 |
| Co-expression Network | WGCNA | Weighted correlation network analysis | Langfelder & Horvath, 2008 |
| Visualization | pheatmap, ggrepel | Heatmaps with clustering, labeled plots | Kolde, 2019; Slowikowski, 2024 |
| Data Manipulation | tidyverse | Data wrangling and transformation | Wickham *et al.*, 2019 |

## Repository Structure

```
rna-seq/
├── README.md                      # This file
├── scripts/
│   ├── 00_install_packages.R      # Dependency installation
│   ├── 01_download_data.R         # GEO data retrieval
│   ├── 02_qc_and_filter.R         # QC plots and gene filtering
│   ├── 03_differential_expression.R # DESeq2 analysis
│   ├── 04_enrichment.R            # gProfiler2 functional enrichment
│   ├── 05_visualization.R         # Publication figures
│   ├── 06_session_summary.R       # Results documentation
│   ├── 07_temporal_patterns.R     # maSigPro temporal clustering
│   └── 08_wgcna_network.R         # WGCNA (exploratory, see Limitations)
├── data/
│   ├── raw/                       # Downloaded data (gitignored)
│   └── processed/                 # RDS files (counts, metadata, VST matrix)
└── results/
    ├── figures/                   # PNG outputs (300 DPI)
    └── tables/                    # CSV result tables
```

## Key Findings

### Differential Expression Analysis

- **Genes analyzed:** 15,893 (after CPM-based filtering)
- **Significant DE genes (30m vs 6m):** 68 total (|log2FC| > 1, padj < 0.05)
  - Upregulated: 66 genes
  - Downregulated: 2 genes
- **PCA:** Clear separation of samples by age group along PC1
- **Top upregulated genes:** C4b, Serpina3n, Gfap (immune/glial activation)

| Comparison | Upregulated | Downregulated | Total |
|------------|-------------|---------------|-------|
| 15m vs 6m  | 5           | 1             | 6     |
| 24m vs 15m | 1           | 1             | 2     |
| 30m vs 24m | 0           | 0             | 0     |
| 30m vs 6m  | 66          | 2             | 68    |

### Functional Enrichment (30m vs 6m upregulated genes)

- **GO Biological Process terms:** 356 significant (FDR < 0.05)
- **KEGG pathways:** 29 significant
- **Top enriched themes:**
  - Immune system processes
  - Complement cascade activation
  - Inflammatory response
  - Glial cell activation

### Temporal Pattern Analysis (maSigPro)

- **Significant genes:** 845 genes with temporal expression patterns
- **Clusters identified:** 5 distinct temporal profiles
  - Cluster 1 (Linear Increase): 222 genes
  - Cluster 2 (Late Response): 148 genes
  - Cluster 3 (Linear Decrease): 132 genes
  - Cluster 4 (Linear Decrease): 181 genes
  - Cluster 5 (Early Response): 162 genes

## Limitations and Methodological Notes

### WGCNA Analysis: Sample Size Considerations

The WGCNA co-expression network analysis script (`08_wgcna_network.R`) was implemented and executed as a demonstration of the methodology. However, **results are not reported** in this analysis due to sample size limitations.

**Rationale:** The dataset contains 16 samples total (n=4 per age group). According to the WGCNA developers' official guidelines:

> *"We do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples... If at all possible, one should have at least 20 samples; as with any analysis methods, more samples usually lead to more robust and refined results."*
>
> — [WGCNA FAQ](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html), Langfelder & Horvath

With only 16 samples (at the absolute minimum threshold), correlation estimates become unstable, and module detection may be driven by sampling noise rather than true biological signal. While the script demonstrates technical proficiency with the WGCNA workflow (soft threshold selection, signed network construction, module-trait correlation, hub gene identification), the resulting modules should be considered **exploratory** and would require validation with a larger cohort.

**Recommendation:** For robust WGCNA results, a minimum of 20-30 samples is advised. The script is retained for methodological reference and can be applied to appropriately powered datasets.

### General Limitations

- Small sample size (n=4 per group) limits statistical power for detecting subtle expression changes
- Bulk RNA-Seq cannot resolve cell-type-specific changes; the original study includes single-cell data for this purpose
- Biological interpretation is intentionally minimal; this analysis focuses on demonstrating the statistical pipeline

## How to Reproduce

### Prerequisites

- R (>= 4.2.0)
- RStudio (recommended) or command-line R
- Internet connection (for GEO data download and gProfiler queries)

### Installation

```bash
# Clone repository
git clone https://github.com/[username]/rna-seq-aging.git
cd rna-seq-aging

# Install R dependencies
Rscript scripts/00_install_packages.R
```

### Run Analysis

Execute scripts sequentially:

```bash
# 1. Download data from GEO
Rscript scripts/01_download_data.R

# 2. Quality control and filtering
Rscript scripts/02_qc_and_filter.R

# 3. Differential expression analysis
Rscript scripts/03_differential_expression.R

# 4. Functional enrichment
Rscript scripts/04_enrichment.R

# 5. Generate publication figures
Rscript scripts/05_visualization.R

# 6. Generate results summary
Rscript scripts/06_session_summary.R

# 7. Temporal pattern analysis
Rscript scripts/07_temporal_patterns.R

# 8. (Optional) WGCNA - exploratory only
Rscript scripts/08_wgcna_network.R
```

### Expected Runtime

- Total pipeline: ~15-30 minutes (depending on hardware)
- GEO download: ~2-5 minutes
- WGCNA (if run): ~5-10 minutes

### Output

Results are written to:
- `results/figures/` - PNG visualizations (300 DPI)
- `results/tables/` - CSV data tables
- `results/session_info.txt` - R session information

## References

### Data Source

Jiang, X., Kang, Y., Pan, L., *et al.* (2024). An atlas of the aging mouse proteome reveals the features of age-related post-transcriptional dysregulation. *Nature Communications*, 15, 8582. https://doi.org/10.1038/s41467-024-52845-x

GEO Accession: [GSE225576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225576)

### Bioinformatics Tools

Conesa, A., Nueda, M. J., Ferrer, A., & Talon, M. (2006). maSigPro: a method to identify significantly differential expression profiles in time-course microarray experiments. *Bioinformatics*, 22(9), 1096-1102. https://doi.org/10.1093/bioinformatics/btl056

Davis, S., & Meltzer, P. S. (2007). GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics*, 23(14), 1846-1847. https://doi.org/10.1093/bioinformatics/btm254

Kolberg, L., Raudvere, U., Kuzmin, I., Vilo, J., & Peterson, H. (2020). gprofiler2 -- an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler. *F1000Research*, 9, 709. https://doi.org/10.12688/f1000research.24956.2

Kolde, R. (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9, 559. https://doi.org/10.1186/1471-2105-9-559

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

Slowikowski, K. (2024). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.9.5. https://CRAN.R-project.org/package=ggrepel

Stephens, M. (2017). False discovery rates: a new deal. *Biostatistics*, 18(2), 275-294. https://doi.org/10.1093/biostatistics/kxw041

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. https://ggplot2.tidyverse.org

Wickham, H., Averick, M., Bryan, J., *et al.* (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686. https://doi.org/10.21105/joss.01686

---

*This analysis was conducted as a portfolio demonstration of RNA-Seq bioinformatics skills.*
