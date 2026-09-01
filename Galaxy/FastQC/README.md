# FastQC

FastQC was used for initial quality control of the raw FASTQ files.

## Purpose

FastQC was performed to evaluate the quality of the sequencing reads and to identify potential sequencing or library preparation issues.

The following quality metrics were inspected:

- Per base sequence quality
- Per sequence quality scores
- Per base sequence content
- Sequence duplication levels
- Adapter content
- Overrepresented sequences

## Results

The 6-month-old samples did not show major adapter contamination and were processed without adapter trimming.

Adapter contamination was detected in the 14-month-old samples. Therefore, adapter trimming and quality filtering were subsequently performed using fastp.

## Input

Paired-end FASTQ files.

## Output

FastQC HTML reports and associated quality-control results.
