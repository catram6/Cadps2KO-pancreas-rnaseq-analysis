# Pancreatic RNA-seq Analysis

## Overview

This repository contains the RNA-seq analysis performed to investigate gene expression changes associated with aging and genotype in mouse pancreatic tissue.

RNA-seq data were obtained from four experimental groups:

- 6-month-old WT
- 6-month-old KO
- 14-month-old WT
- 14-month-old KO

The analysis was performed using Galaxy and R.

## Experimental design

Four pairwise comparisons were performed:

1. 6-month WT vs KO
2. 14-month WT vs KO
3. WT 6-month vs 14-month
4. KO 6-month vs 14-month

## RNA-seq analysis workflow

```text
FASTQ files
↓
FastQC
↓
fastp (14-month samples)
↓
RNA STAR
↓
MarkDuplicates
↓
featureCounts
↓
Count matrix generation
↓
DESeq2
↓
Differentially expressed genes (DEGs)
↓
Volcano plots
↓
DEG overlap analysis
↓
7 gene subgroups
↓
Metascape enrichment analysis
↓
ChIP-Atlas enrichment analysis
```

### Quality control and preprocessing

FASTQ files were initially evaluated using FastQC.

For the 14-month-old samples, adapter contamination was detected by FastQC. Therefore, adapter trimming and quality filtering were performed using fastp.

The following parameters were used for paired-end reads:

```bash
fastp -i input_R1.fq.gz -I input_R2.fq.gz -o trimmed_R1.fq.gz -O trimmed_R2.fq.gz --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 36 -w 4
```

The 6-month-old samples were processed without fastp because no major adapter contamination was detected.

### Alignment and read counting

Reads were aligned to the mouse reference genome using RNA STAR.

PCR duplicates were removed using MarkDuplicates.

Gene-level read counts were obtained using featureCounts.

Count files were combined using Column Join in Galaxy to generate count matrices for the four differential expression comparisons.

## Differential expression analysis

Differential expression analysis was performed using DESeq2.

The following workflow was applied to each of the four comparisons:

```r
library(DESeq2)

samples <- read.csv("sample_table.csv")

samples$condition <- factor(samples$condition, levels = c("reference", "comparison"))

counts_df <- read.csv("count_matrix.csv", row.names = 1)

count_matrix <- as.matrix(counts_df)

storage.mode(count_matrix) <- "numeric"

count_matrix <- na.omit(count_matrix)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = samples, design = ~ condition)

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds)
```

### Count matrix preprocessing

The original count matrices contained an additional non-numeric row.

During conversion of the count matrix to numeric format, this resulted in NA values.

These rows were removed using:

```r
count_matrix <- na.omit(count_matrix)
```

This preprocessing step was applied to the count matrices used for the DESeq2 analyses.

### DEG identification

DEGs were identified based on adjusted p-value and log2 fold change.

Upregulated genes were defined as:

```r
padj < 0.05 & log2FoldChange > threshold
```

Downregulated genes were defined as:

```r
padj < 0.05 & log2FoldChange < -threshold
```

Gene IDs were added using:

```r
res$gene <- rownames(res)
up$gene <- rownames(up)
down$gene <- rownames(down)
```

The thresholds used for each comparison were:

| Comparison | Adjusted p-value | log2FC threshold |
|---|---:|---:|
| 6-month WT vs KO | < 0.05 | ｜log2FC｜ > 1 |
| 14-month WT vs KO | < 0.05 | ｜log2FC｜ > 1 |
| WT 6-month vs 14-month | < 0.05 | ｜log2FC｜ > 2 |
| KO 6-month vs 14-month | < 0.05 | ｜log2FC｜ > 2 |

The log2 fold-change thresholds were increased for some comparisons because the number of DEGs exceeded 3,000, which prevented downstream analysis using Metascape.

### DEG results

The numbers of identified DEGs were:

| Comparison | Upregulated | Downregulated | Total |
|---|---:|---:|---:|
| 6-month WT vs KO | 1,871 | 1,813 | 3,684 |
| 14-month WT vs KO | 701 | 448 | 1,149 |
| WT 6-month vs 14-month | 1,490 | 1,090 | 2,580 |
| KO 6-month vs 14-month | 1,221 | 835 | 2,056 |

The complete DESeq2 results and DEG lists are stored in the results/ directory.

## Volcano plot

Volcano plots were generated to visualize differential gene expression.

The plots display:

- log2 fold change on the x-axis
- statistical significance on the y-axis
- significantly upregulated genes
- significantly downregulated genes

The volcano plots were generated from the DESeq2 results.

## DEG overlap analysis

The DEG lists obtained from the four comparisons were compared to identify genes shared between different experimental comparisons.

Genes were classified according to their occurrence patterns across the four DEG datasets.

Based on these patterns, genes were divided into seven subgroups.

## Functional enrichment analysis

Functional enrichment analysis was performed for each of the seven gene subgroups.

### Metascape

Each gene subgroup was analyzed using Metascape to investigate:

- enriched biological processes
- pathways
- functional categories
- Gene Ontology terms

### ChIP-Atlas

ChIP-Atlas enrichment analysis was performed for each gene subgroup to investigate potential transcription factors and regulatory mechanisms associated with the genes.

The results from Metascape and ChIP-Atlas were used to investigate the biological characteristics and potential regulatory mechanisms of each gene subgroup.

## Repository structure

```text
.
├── README.md
├── R/
├── galaxy/
├── results/
└── docs/
```

## Analysis status

The computational analysis described above has been completed.

Biological interpretation and experimental validation are ongoing.
