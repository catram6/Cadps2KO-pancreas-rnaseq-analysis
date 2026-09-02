# Analysis Overview

## Research objective

This study aimed to investigate gene expression changes associated with pancreatic aging and the effect of Cadps2 knockout (KO) on age-associated pancreatic changes.

RNA-seq data from mouse pancreatic tissue were analyzed to identify differentially expressed genes (DEGs) associated with genotype, aging, and their interaction.

## Experimental design

Four experimental groups were analyzed:

- 6-month-old WT
- 6-month-old KO
- 14-month-old WT
- 14-month-old KO

These groups were used to perform four pairwise differential expression analyses.

## Differential expression comparisons

The following four comparisons were performed:

| Comparison | Biological question |
|---|---|
| 6-month WT vs KO | Effect of Cadps2 KO at 6 months |
| 14-month WT vs KO | Effect of Cadps2 KO at 14 months |
| WT 6-month vs 14-month | Age-associated changes in WT |
| KO 6-month vs 14-month | Age-associated changes in KO |

## RNA-seq analysis

Raw paired-end FASTQ files were processed using the following workflow:

```text
FASTQ
↓
FastQC
↓
fastp
↓
RNA STAR
↓
MarkDuplicates
↓
featureCounts
↓
Count matrix
↓
DESeq2
```

Adapter contamination was detected in the 14-month-old samples by FastQC. Therefore, fastp was used for adapter trimming and quality filtering of these samples.

The 6-month-old samples did not undergo fastp preprocessing.

Reads were aligned to the mouse mm10 reference genome using RNA STAR with the GENCODE vM10 annotation.

PCR duplicates were removed using MarkDuplicates.

Gene-level read counts were generated using featureCounts.

## Differential expression analysis

Differential expression analysis was performed using DESeq2.

Genes were considered DEGs when they satisfied both an adjusted p-value threshold and a log2 fold-change threshold.

The thresholds used were:

| Comparison | padj | ｜log2FC｜ |
|---|---:|---:|
| 6-month WT vs KO | < 0.05 | > 2 |
| 14-month WT vs KO | < 0.05 | > 2 |
| WT 6-month vs 14-month | < 0.05 | > 2 |
| KO 6-month vs 14-month | < 0.05 | > 2 |

The same log2 fold-change threshold of 2 was ultimately used for all four comparisons.

## DEG overlap analysis

The DEG lists from the four comparisons were integrated to identify genes shared between different biological conditions.

For the overlap analysis, genes were classified according to the combinations of comparisons in which they were identified as DEGs.

For comparisons involving upregulated and downregulated genes, the direction of change was taken into account so that genes with consistent directional changes were grouped together.

No subgroup corresponding to genes shared only between the KO 6-month vs 14-month comparison and the 6-month WT vs KO comparison was identified.

Based on the overlap patterns, the DEGs were classified into seven subgroups:

1. KO and WT aging-associated DEGs
2. DEGs shared by KO aging, WT aging, and 14-month WT vs KO
3. DEGs shared by KO aging, WT aging, and 6-month WT vs KO
4. DEGs shared by all four comparisons
5. DEGs shared by KO aging, 14-month WT vs KO, and WT aging
6. DEGs shared by KO aging and 14-month WT vs KO
7. DEGs identified only in KO aging

The seven subgroups were subsequently used for functional and regulatory enrichment analyses.

## Functional enrichment analysis

Each DEG subgroup was analyzed separately using Metascape and ChIP-Atlas.

### Metascape

Metascape was used to investigate functional characteristics of each subgroup.

### ChIP-Atlas

ChIP-Atlas enrichment analysis was performed to investigate potential transcriptional regulators and regulatory mechanisms associated with each DEG subgroup.

## Volcano plot analysis

Volcano plots were generated to visualize the differential expression results.

For the KO aging comparison, a customized volcano plot was also generated to highlight genes associated with pancreatic phenotypes, including:

- Acinar / secretory granule
- Exocytosis machinery
- Fibrosis / immune infiltration
- ADM / ductal

Cadps2 was specifically excluded from this visualization when highlighting these phenotype-associated genes.

## Current status

The analysis of the current datasets has been completed.

The overall research project is ongoing, with biological interpretation of the DEG subgroups and enrichment analysis results in progress. Additional datasets will be incorporated, and experimental validation will be performed in future studies.
