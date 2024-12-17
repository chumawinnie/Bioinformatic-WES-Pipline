# bioinformatic-WES-pipline
Whole-Exon-Sequncing 🔬🧬📊💻🧪📈✨
# Bioinformatics WES Pipeline

## Overview
This pipeline processes Whole-Exon Sequencing (WES) data to generate variant calls and annotations.

## Requirements
- Python 3.8+
- Perl
- Required tools: `bwa`, `samtools`, `vcf2maf`

## Pipeline Steps
1. **Pre-processing**:
   - Align reads using BWA.
   - Sort and index BAM files.
2. **Variant Calling**:
   - Use Mutect2 for somatic variant detection.
3. **Annotation**:
   - Annotate VCF files with VEP and convert to MAF format.

## Usage
```bash
bash merged-sample-WES.sh
bash single-sample-WES.sh
