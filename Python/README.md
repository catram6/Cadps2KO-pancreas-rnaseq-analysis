# Python Analysis

This directory contains Python scripts used for visualization of differential expression analysis results.

## Volcano plots

Python was used to generate Volcano plots from the DESeq2 results.

The plots visualize:

- log2 fold change on the x-axis
- -log10 adjusted p-value on the y-axis
- significantly upregulated genes
- significantly downregulated genes

The main libraries used were:

- pandas
- NumPy
- Matplotlib
- openpyxl
- JSON

## Standard Volcano plots

Volcano plots were generated for three comparisons:

1. 6-month-old KO vs WT
2. 14-month-old KO vs WT
3. WT 14-month vs 6-month

The DEG criteria used for visualization correspond to the thresholds used in the DESeq2 analysis.

### DEG thresholds

For genotype comparisons:

`padj < 0.05 & |log2FC| > 1`

For aging comparisons:

`padj < 0.05 & |log2FC| > 2`

## KO aging Volcano plot

A customized Volcano plot was generated for the KO aging comparison.

The plot highlights genes related to pancreatic phenotypes and biological processes.

The highlighted genes were classified into the following categories:

- Acinar / secretory granule
- Exocytosis machinery
- Fibrosis / immune infiltration
- ADM / ductal

Cadps2 was excluded from the highlighted category plot to allow visualization of other phenotype-associated genes.

## Output

The generated Volcano plots are stored in the `results/` directory.

The Python scripts used to generate the plots are stored in this directory.

## Reproducibility

The scripts contain the parameters used to generate the figures, including:

- significance thresholds
- log2 fold-change thresholds
- plotting limits
- figure dimensions
- output formats
- highlighted gene categories

The scripts can be used to reproduce the Volcano plots from the corresponding DESeq2 result files.
