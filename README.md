

---

### Full Updated `README.md`

Copy and paste this into your `README.md` file:

```markdown
# Bioinformatics WES Pipeline

## Overview
This pipeline processes Whole-Exome Sequencing (WES) data to perform variant discovery, copy number variation (CNV) analysis, single nucleotide polymorphism (SNP) detection, and variant annotation. It integrates several bioinformatics tools to analyze tumor-normal paired samples for downstream interpretation.

The pipeline was designed to support:
- **Somatic Variant Detection** (Mutect2)
- **Copy Number Variations (CNV)** (CNVkit and Control-FREEC)
- **Loss of Heterozygosity (LOH)** (Control-FREEC)
- **Single Nucleotide Polymorphisms (SNP)** analysis
- **Variant Annotation** (VEP, vcf2maf)
- **Functional Enrichment** (Gene Ontology, KEGG pathway)

This pipeline is structured to be reproducible and efficient for cancer genomics studies.

---

## Requirements
The pipeline requires the following tools and dependencies:

### Programming Languages
- **Python 3.8+**
- **Perl 5+**

### Bioinformatics Tools
- **BWA**: Align reads to the reference genome.
- **Samtools**: Process and sort BAM files.
- **Mutect2**: Somatic variant calling.
- **CNVkit**: Copy Number Variation analysis.
- **Control-FREEC**: CNV and LOH detection.
- **VEP**: Variant Effect Predictor for annotation.
- **vcf2maf**: Convert VCF files to MAF format.

### Other Requirements
- Reference genome (e.g., `hg19` or `GRCh37`).
- BED files for exome target intervals.

---

## Pipeline Steps

### 1. **Pre-processing**
- Align paired-end reads to the reference genome using **BWA**.
- Sort and index the resulting BAM files with **Samtools**.
- Mark duplicates and recalibrate base quality scores (optional).

### 2. **Somatic Variant Calling**
- Use **Mutect2** (GATK) to identify somatic single nucleotide variants (SNVs) and indels between tumor and normal samples.
- Generate VCF files with filtered variants.

### 3. **Copy Number Variation (CNV) Analysis**
- Use **CNVkit** to analyze CNVs:
  - Detect gains and losses of genomic regions.
  - Generate plots (scatter, diagram).
- Run **Control-FREEC** for CNV and **Loss of Heterozygosity (LOH)** detection:
  - Tumor-normal comparison for regions of allelic imbalance.

### 4. **SNP Analysis**
- Extract SNPs from the Mutect2 VCF files for detailed analysis.

### 5. **Variant Annotation**
- Annotate filtered variants using **VEP (Variant Effect Predictor)**.
- Convert VCF files to MAF format for downstream analyses using **vcf2maf**.

### 6. **Post-Processing and Functional Enrichment**
- Perform **Gene Ontology (GO)** and **KEGG pathway** enrichment analysis on the identified variants.

---

## Usage

### Step 1: Run Pre-processing and Variant Calling
```bash
bash merged-sample-WES.sh
bash single-sample-WES.sh
```

### Step 2: Copy Number Variation (CNV) Analysis
For CNVkit:
```bash
cnvkit.py batch Tumor.bam -n Normal.bam -t targets.bed --scatter --diagram
```

For Control-FREEC:
- Configure your Control-FREEC `.conf` file.
- Run Control-FREEC:
```bash
freec -conf controlfreec.conf
```

### Step 3: Variant Annotation
Annotate VCF files with VEP:
```bash
vep --input_file Tumor_vs_Normal.vcf --output_file annotated_output.vcf --species homo_sapiens --assembly GRCh37
```

Convert VCF to MAF:
```bash
perl vcf2maf.pl --input-vcf annotated_output.vcf --output-maf output.maf --ref-fasta hg19.fa
```

---

## Outputs

1. **Variant Calling**:
   - Filtered somatic variants (VCF).
   - Summary statistics.

2. **CNV Analysis**:
   - CNV plots (`.png`, `.pdf`).
   - CNV calls (`.cns`, `.cnr` files).

3. **LOH Analysis**:
   - Regions of allelic imbalance.

4. **Annotation**:
   - Annotated VCF files.
   - MAF files for downstream analysis.

5. **Reports**:
   - Gene lists for functional enrichment (GO, KEGG).

---

## Workflow Diagram

![Pipeline Workflow](link_to_diagram.png)

---

## Future Plans
- Integrate Tumor Mutational Burden (TMB) calculation.
- Add tools for **mutation signature analysis** (e.g., SigProfilerExtractor).

---

## References
- GATK: https://gatk.broadinstitute.org
- CNVkit: https://cnvkit.readthedocs.io
- Control-FREEC: http://boevalab.com/FREEC/
- VEP: https://www.ensembl.org/info/docs/tools/vep/index.html

---

## Contact
For questions or issues, feel free to contact me at: **chumawinnie@gmail.com**
```

---

### What This Includes:
1. **All your pipeline steps**: Preprocessing, variant calling, CNVs, LOH, SNPs, and annotation.
2. **Usage instructions**: Commands for running the scripts and tools.
3. **Outputs**: Clear description of the results.
4. **Workflow diagram** placeholder (you can add a diagram later).
5. **Future plans**: Highlights TMB and mutation signature analysis.

---

### What to Do Next:
1. Paste this content into your `README.md` file.
2. Save and push the changes to GitHub:
   ```bash
   git add README.md
   git commit -m "Update README with full pipeline description"
   git push origin master
   ```




