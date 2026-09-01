# DEG Subgroup Analysis

## Purpose

The DEG lists obtained from the four differential expression analyses were integrated to identify genes shared between different biological conditions.

The four DEG datasets were:

1. 6-month WT vs KO
2. 14-month WT vs KO
3. WT 6-month vs 14-month
4. KO 6-month vs 14-month

Genes were classified according to the combinations of comparisons in which they were identified as DEGs.

For genes that were identified as both upregulated and downregulated in different comparisons, the direction of differential expression was considered when determining whether genes were shared.

## DEG direction

The direction of differential expression was defined according to the following comparisons:

- 6-month WT vs KO
- 14-month WT vs KO
- WT 6-month vs 14-month
- KO 6-month vs 14-month

For the aging comparisons, the comparison represents 14-month-old samples relative to 6-month-old samples.

For the genotype comparisons, the comparison represents KO relative to WT.

Genes were grouped only when their DEG status and direction were consistent with the subgroup definition.

## Seven DEG subgroups

The overlap analysis resulted in seven gene subgroups.

| Subgroup | Shared DEG comparisons | Total genes | UP in KO aging | DOWN in KO aging |
|---|---|---:|---:|---:|
| Subgroup 1 | KO aging + WT aging | 757 | 407 | 350 |
| Subgroup 2 | KO aging + WT aging + 14-month WT vs KO | 28 | 27 | 1 |
| Subgroup 3 | KO aging + WT aging + 6-month WT vs KO | 162 | 54 | 108 |
| Subgroup 4 | KO aging + WT aging + 14-month WT vs KO + 6-month WT vs KO | 31 | 17 | 14 |
| Subgroup 5 | KO aging + 14-month WT vs KO + 6-month WT vs KO | 3 | 2 | 1 |
| Subgroup 6 | KO aging + 14-month WT vs KO | 161 | 152 | 9 |
| Subgroup 7 | KO aging only | 914 | 562 | 352 |

## Subgroup descriptions

### Subgroup 1: Aging-associated DEGs shared by WT and KO

Genes that showed differential expression during aging in both WT and KO mice.

These genes represent aging-associated transcriptional changes that occur in both genotypes.

### Subgroup 2: Aging-associated DEGs shared with the 14-month genotype comparison

Genes identified as aging-associated DEGs in both WT and KO and also differentially expressed between WT and KO at 14 months.

These genes may represent aging-associated changes that are particularly associated with genotype differences in aged animals.

### Subgroup 3: Aging-associated DEGs shared with the 6-month genotype comparison

Genes identified as aging-associated DEGs in both WT and KO and also differentially expressed between WT and KO at 6 months.

These genes may represent genes whose genotype-associated differences are already detectable at 6 months and are also altered during aging.

### Subgroup 4: DEGs shared by all four comparisons

Genes identified as DEGs in all four differential expression comparisons.

These genes showed differential expression associated with both genotype and aging across the experimental conditions.

### Subgroup 5: DEGs shared by KO aging and both genotype comparisons

Genes identified as DEGs in KO aging and in both the 6-month and 14-month WT vs KO comparisons, but not in WT aging.

These genes may represent genotype-associated changes that are particularly pronounced in the KO background.

### Subgroup 6: DEGs shared by KO aging and the 14-month genotype comparison

Genes identified as DEGs in KO aging and in the 14-month WT vs KO comparison, but not in the other comparisons.

These genes may represent changes associated with aging in KO mice that are also associated with genotype differences specifically at 14 months.

### Subgroup 7: KO-aging-specific DEGs

Genes identified as DEGs only in the KO 6-month vs 14-month comparison.

These genes represent transcriptional changes associated with aging that were detected specifically in the KO background.

## Absence of a KO-aging + 6-month genotype-only subgroup

No DEG subgroup consisting exclusively of genes shared between:

- KO 6-month vs 14-month
- 6-month WT vs KO

was identified.

Therefore, this overlap pattern was not included as a separate subgroup.

## Downstream analysis

Each of the seven DEG subgroups was analyzed independently using:

- Metascape
- ChIP-Atlas

These analyses were performed to investigate the functional characteristics and potential regulatory mechanisms associated with each subgroup.
