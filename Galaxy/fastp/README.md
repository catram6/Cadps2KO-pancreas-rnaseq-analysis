# fastp

Adapter trimming and quality filtering were performed using fastp for the 14-month-old RNA-seq samples.

## Reason for trimming

FastQC identified adapter contamination in the 14-month-old samples.

Therefore, adapter trimming and quality filtering were performed before alignment.

The 6-month-old samples were not processed using fastp because no major adapter contamination was detected.

## Tool

fastp

## Data type

Paired-end RNA-seq reads.

## Parameters

The following parameters were used:

```bash
--detect_adapter_for_pe
--qualified_quality_phred 20
--length_required 36
-w 4
```

## Command

The following command was used for the analysis:

```bash
fastp -i input_R1.fq.gz -I input_R2.fq.gz -o trimmed_R1.fq.gz -O trimmed_R2.fq.gz --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 36 -w 4
```

## Parameter description

- `--detect_adapter_for_pe`: Automatically detects adapters in paired-end reads.
- `--qualified_quality_phred 20`: Requires a base quality score of at least Phred 20.
- `--length_required 36`: Removes reads shorter than 36 bases after filtering.
- `-w 4`: Uses 4 processing threads.

## Output

fastp generates:

- Trimmed paired-end FASTQ files
- HTML quality-control report
- JSON quality-control report

The HTML and JSON reports were generated during processing.
