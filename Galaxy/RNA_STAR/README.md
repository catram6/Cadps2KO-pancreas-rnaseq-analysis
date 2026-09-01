# RNA STAR

RNA STAR was used to align RNA-seq reads to the mouse reference genome.

## Purpose

RNA STAR was used to perform spliced alignment of RNA-seq reads to the mouse reference genome and generate BAM files for downstream analysis.

## Input

Paired-end RNA-seq FASTQ files.

For the 14-month-old samples, FASTQ files after fastp preprocessing were used.

For the 6-month-old samples, the original FASTQ files were used.

## Parameters

The following settings were used:

- Read type: Paired-end
- Reference genome: mm10
- Gene model: gencode.vM10.annotation.gtf
- Length of the genomic sequence around annotated junctions: 149

Other parameters were left at their default settings.

## Alignment

RNA STAR was used to perform spliced alignment of the RNA-seq reads.

The resulting BAM files were used for downstream duplicate removal and gene-level read counting.

## Output

BAM files containing the aligned RNA-seq reads.

These BAM files were subsequently processed using MarkDuplicates and featureCounts.

## Workflow

```text
FASTQ
↓
RNA STAR
↓
BAM
↓
MarkDuplicates
↓
featureCounts
```
