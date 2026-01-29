# GSE225576 Mouse Brain Aging RNA-seq - Results Index

## Project Overview

This analysis examines differential gene expression in mouse brain tissue across four age groups:
- 6 months (young adult)
- 15 months (middle age)
- 24 months (early aged)
- 30 months (aged)

## Key Findings

- **Total samples analyzed:** 16
- **Genes after filtering:** 15,893
- **Differential expression (30m vs 6m):** 66 up, 2 down
- **Temporal patterns (maSigPro):** 845 significant genes, 5 clusters
- **Major biological themes:**
  - Immune response and inflammation (strongest enrichment)
  - Complement system activation
  - Glial cell activation (astrogliosis)
  - Synaptic transmission downregulation

## Generated Figures

### QC and DE Analysis
- `DE_summary.png`: Summary of differential expression across all comparisons
- `enrichment_dotplot.png`: GO and KEGG enrichment analysis results
- `heatmap_top50_annotated.png`: Heatmap with gene functional cluster annotations
- `heatmap_top50.png`: Heatmap of top 50 differentially expressed genes
- `marker_genes_significance.png`: Marker genes with statistical significance brackets
- `marker_genes.png`: Expression trends for aging marker genes
- `qc_correlation.png`: Sample correlation heatmap
- `qc_library_size.png`: Library size distribution across samples
- `qc_pca.png`: PCA plot showing sample clustering by age group
- `volcano_30m_vs_6m_polished.png`: Volcano plot with ggrepel labels (polished)
- `volcano_30m_vs_6m.png`: Volcano plot: 30m vs 6m differential expression

### Temporal Pattern Analysis (maSigPro)
- `maSigPro_clusters.png`: Temporal expression pattern clusters (5 clusters)
- `maSigPro_heatmap.png`: Heatmap of significant temporal genes

## Generated Tables

### Differential Expression
- `analysis_summary.csv`: Comprehensive analysis summary metrics
- `DE_15m_vs_6m.csv`: Differential expression results: 15m vs 6m
- `DE_24m_vs_15m.csv`: Differential expression results: 24m vs 15m
- `DE_30m_vs_24m.csv`: Differential expression results: 30m vs 24m
- `DE_30m_vs_6m.csv`: Differential expression results: 30m vs 6m

### Functional Enrichment
- `gProfiler_GO_KEGG.csv`: Functional enrichment results (GO + KEGG)

### Temporal Analysis (maSigPro)
- `maSigPro_significant_genes.csv`: Genes with significant temporal patterns
- `maSigPro_clusters.csv`: Cluster assignments for significant genes
- `cluster_profiles.csv`: Mean expression profiles for each cluster

## Methodology

1. **Quality Control:** Library size, PCA, and correlation analysis
2. **Filtering:** Genes with CPM > 1 in at least 3 samples retained
3. **Normalization:** DESeq2 variance stabilizing transformation (VST)
4. **Differential Expression:** DESeq2 with ashr LFC shrinkage
5. **Enrichment:** gProfiler2 (GO:BP + KEGG, FDR corrected)
6. **Temporal Patterns:** maSigPro polynomial regression (degree=2) with k-means clustering

## Note on WGCNA

WGCNA co-expression network analysis was implemented (see `scripts/08_wgcna_network.R`) but results are not reported due to sample size limitations (n=16). The WGCNA developers recommend a minimum of 20 samples for robust module detection.

## Session Information

See `session_info.txt` for detailed R environment information.
