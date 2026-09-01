# R Analysis

This directory contains the R scripts used for differential expression analysis of the pancreatic RNA-seq data.

## Differential expression analysis

Differential expression analysis was performed using DESeq2 for four pairwise comparisons:

1. 6-month-old WT vs KO
2. 14-month-old WT vs KO
3. WT 6-month-old vs 14-month-old
4. KO 6-month-old vs 14-month-old

The four scripts correspond to these comparisons:

- `01_DESeq2_6m_WT_vs_KO.R`
- `02_DESeq2_14m_WT_vs_KO.R`
- `03_DESeq2_WT_aging.R`
- `04_DESeq2_KO_aging.R`

## DESeq2 workflow

Each analysis followed the general workflow:

```text
Count matrix  
↓  
Read count data and sample information  
↓  
Convert count matrix to numeric format  
↓  
Remove non-numeric rows  
↓  
Create DESeqDataSet  
↓  
Filter low-count genes  
↓  
DESeq2 normalization and dispersion estimation  
↓  
Differential expression testing  
↓  
DEG extraction  
↓  
Export results as CSV
```

## Count matrix preprocessing

The count matrices generated from the Galaxy analysis contained an additional non-numeric row.

The count matrices were converted to numeric format using:

```r
storage.mode(count_matrix) <- "numeric"
```

This generated NA values for the non-numeric row.

The resulting NA-containing rows were removed using:

```r
count_matrix <- na.omit(count_matrix)
```

Genes with a total count of less than 10 across all samples were then excluded:

```r
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]
```

## DEG criteria

DEGs were identified using an adjusted p-value threshold of 0.05.

The log2 fold-change threshold differed between genotype comparisons and aging comparisons.

### Genotype comparisons: WT vs KO

For the 6-month and 14-month WT vs KO comparisons, DEGs were defined as:

**Upregulated:**

`padj < 0.05 & log2FoldChange > 1`

**Downregulated:**

`padj < 0.05 & log2FoldChange < -1`

Therefore, the DEG criterion was:

`|log2FoldChange| > 1`

### Aging comparisons: 6-month vs 14-month

For the WT and KO aging comparisons, DEGs were defined using a more stringent fold-change threshold:

**Upregulated:**

`padj < 0.05 & log2FoldChange > 2`

**Downregulated:**

`padj < 0.05 & log2FoldChange < -2`

Therefore, the DEG criterion was:

`|log2FoldChange| > 2`

## DEG criteria summary

| Comparison | Adjusted p-value | log2FC threshold |
|---|---:|---:|
| 6-month WT vs KO | < 0.05 | ｜log2FC｜ > 1 |
| 14-month WT vs KO | < 0.05 | ｜log2FC｜ > 1 |
| WT 6-month vs 14-month | < 0.05 | ｜log2FC｜ > 2 |
| KO 6-month vs 14-month | < 0.05 | ｜log2FC｜ > 2 |

## Output

The DESeq2 scripts generate:

- complete DESeq2 results
- upregulated DEG lists
- downregulated DEG lists

The resulting files are stored in the `results/` directory.

## Downstream analysis

The DESeq2 results were subsequently used for:

- Volcano plot visualization
- DEG overlap analysis
- Seven DEG subgroup classification
- Metascape enrichment analysis
- ChIP-Atlas enrichment analysis

## Notes

The DESeq2 analysis was performed as part of an ongoing study of pancreatic aging and disease-associated gene expression changes.

The biological interpretation and experimental validation of the results are ongoing.
