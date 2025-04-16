
---

# **Bioinformatics WES Pipeline** 💉🩸🧪🔬🧬💻📈📊✨  

## **Overview**  
This pipeline processes **Whole-Exome Sequencing (WES)** data for comprehensive cancer genomics analysis. It covers variant discovery, copy number variation (CNV), single nucleotide polymorphisms (SNPs), indel detection, and variant annotation. Additionally, the pipeline supports **complex biomarkers analysis** 🧪🧬✨💊🔍 such as:  

- **Homologous Recombination Deficiency (HRD)**  
- **HLA Typing** for immune system markers  
- **Epitope Prediction** for neoantigens
- **Microsatellite Instability (MSI)**  
- **Mutational Signature Detection**
- **Tumour Mutational Burden**

The pipeline integrates widely-used bioinformatics tools and supports tumor-normal paired samples for downstream interpretation, aiding molecular tumor boards and cancer research.

---

## **Key Features**  ✨🔑📋💡

- **Somatic Variant Detection** (Mutect2)  
- **Copy Number Variations (CNV)**: CNVkit and Control-FREEC  
- **Loss of Heterozygosity (LOH)**: Control-FREEC  
- **Single Nucleotide Polymorphism (SNP)** detection  
- **Cellularity and Ploidy Estimation**: Sequenza  
- **HRD Analysis**: LOH, TAI, LST scores (Sequenza, scarHRD)  
- **Variant Annotation**: VEP and VCF-to-MAF conversion
- **Microsatellite Instability (MSI)**: MSIsensor-pro 
- **Tumor Mutational Burden** : python custom script 
- **Mutational Signature Analysis**: SigProfilerExtractor  
- **Human Leukocyte Antigen (HLA) Typing**: OptiType  
- **Epitope Prediction**: pVACseq for neoantigen discovery  
- **Functional Enrichment Analysis**: Gene Ontology (GO), KEGG Pathways  

---

## **Requirements** 💻💻👨‍💻👩‍💻📟  

### **Programming Languages**  
- 🐍 (Python) 3.8+  
- 🐪 Perl 5+  
- 📊R 4.4.2  
- 🐚Bash Shell 5.2  
- 🚀Nextflow 24.10.2

### **Bioinformatics Tools** 💻👨‍💻👩‍💻📟 🛠️🔧🪛⚙️
1. **Preprocessing & Alignment**:
   - Fastqc: Quality-control Processing 
   - BWA: Align reads to the reference genome.  
   - Samtools: Process and sort BAM files.  

3. **Variant Calling & Annotation**:  
   - Mutect2: Detect somatic SNVs and indels.  
   - VEP (Variant Effect Predictor): Annotate variants.  
   - vcf2maf: Convert VCF files to MAF format.  

4. **Copy Number Variation (CNV)**:  
   - CNVkit: CNV detection, scatter plots, and genome diagrams.  
   - Control-FREEC: CNV and LOH analysis.  

5. **HRD Analysis**:  
   - Sequenza: Estimate tumor cellularity, ploidy, and HRD metrics.  
   - scarHRD: Calculate LOH, TAI, and LST scores for HRD detection.  

6. **Immune System Analysis**:  
   - **OptiType**: HLA Typing to identify Class I HLA alleles.  
   - **pVACseq**: Predict epitopes (neoantigens) for MHC Class I & II.  

7. **Mutational Signature Analysis**:  
   - **SigProfilerExtractor**: Extract and analyze mutational signatures.
     
8. - **Tumor Mutational Burden**: calculate / measures the total number of somatic mutations per megabase (Mb) 

### **Other Requirements**  
- Reference genome: hg19 or GRCh37  
- BED files: Target capture regions for exome sequencing  

---

## **Pipeline Steps**  🪜➡️🔄📋

### **1. Pre-processing**  
- Align paired-end reads to the reference genome using **BWA**.  
- Sort and index BAM files using **Samtools**.  
- Mark duplicates and recalibrate quality scores using **Picard**.  

### **2. Somatic Variant Calling**  
- Detect **somatic SNVs** and **indels** using **Mutect2**.  
- Filter low-confidence variants and generate a clean **VCF file**.  

### **3. Copy Number Variation (CNV) Analysis**  
- **CNVkit**:  
  - Detect gains and losses of genomic regions.  
  - Generate scatter plots and diagrams for CNVs.  
- **Control-FREEC**:  
  - Perform CNV detection and **LOH** analysis.  
  - Identify allelic imbalances in tumor-normal samples.

  ### **4. MSIsensor-pro(MSI)**
- **MSIsensor-pro**:
  - Detect Microsatellite Instability (MSI)
  - Determines whether a tumor is microsatellite stable (MSS)

### **5. Tumor Mutational Burden**
- **Tumor Mutational Burden**:
  - calculate total number of mutation per mergabase (tmb)
  - Determines whether a tumor is high or low 

### **6. HRD Analysis**  
- **Sequenza**: Estimate tumor cellularity and ploidy.  
- **scarHRD**:  
   - Calculate metrics like **LOH (Loss of Heterozygosity)**, **TAI (Telomeric Allelic Imbalance)**, and **LST (Large-scale Transitions)**.  
   - These metrics quantify **Homologous Recombination Deficiency (HRD)**, a key biomarker in DNA repair-deficient tumors.  

### **7. HLA Typing**  
- **OptiType**: Analyze immune system HLA alleles from sequencing data to identify Class I HLA genes (e.g., HLA-A, HLA-B).  

### **8. Epitope Prediction**  
- **pVACseq**:  
   - Predict neoantigens (epitopes) that bind to **MHC Class I & II** molecules.  
   - Helps identify targets for cancer vaccines or immunotherapies.  

### **9. Mutational Signature Analysis**  
- **SigProfilerExtractor**: Analyze mutational patterns to identify **mutation signatures** associated with processes like:  
  - Aging  
  - Smoking  
  - DNA repair defects  

### **10. Variant Annotation**  
- Annotate variants using **VEP**.  
- Convert annotated VCF files to MAF format using **vcf2maf**.  

### **11. Functional Enrichment**  
- Perform **Gene Ontology (GO)** and **KEGG pathway** analysis on annotated gene lists.  

---

## **Usage**  📖⚙️🛠️💡

### **Step 1: Pre-processing and Variant Calling**  
```bash
bash merged-sample-WES.sh  
bash single-sample-WES.sh  
```  

### **Step 2: CNV Analysis**  

**CNVkit**:  
```bash
cnvkit.py batch Tumor.bam -n Normal.bam -t targets.bed --scatter --diagram  
```
**Control-FREEC**:  
```bash
freec -conf controlfreec.conf  
```   

### **Step 3: **MSIsensor-pro**
```bash
msisensor-pro msi \
  -d reference.microsatellites.bed \
  -t Tumor.bam \
  -n Normal.bam \
  -o output_prefix \
  -b bed_file_of_interest.bed \
  -c 10
```

 

### **Step 4: HRD Analysis**  
**Sequenza**:  
```bash
Rscript sequenza_analysis.R  
```  

### **Step 5: HLA Typing and Epitope Prediction**  
```bash
python OptiTypePipeline.py  
pvacseq run input.vcf output_dir hg19  
```  

### **Step 6: Variant Annotation**  
**VEP**:  
```bash
vep --input_file input.vcf --output_file annotated_output.vcf --species homo_sapiens --assembly GRCh37  
```  
**VCF to MAF**:  
```bash
perl vcf2maf.pl --input-vcf annotated_output.vcf --output-maf output.maf --ref-fasta hg19.fa  
```  

---

## **Outputs**  📖💡📈📊✨

1. **Variant Calling**: Filtered VCF files and stats  
2. **CNV Analysis**: CNV plots (.png, .pdf) and call files  
3. **LOH Analysis**: Regions of allelic imbalance
4. **MSIsensor-pro**: Detect Microsatellite Instability (MSI)
5. **HRD Metrics**: LOH, TAI, and LST scores  
6. **HLA Typing**: HLA alleles list  
7. **Epitope Prediction**: Neoantigen lists  
8. **Mutational Signatures**: Signature plots and matrices  
9. **Annotation**: Annotated VCF and MAF files  

---

## **Future Plans**  📅🔮➡️📈📊✨
- **Tumor Mutational Burden (TMB)** calculation  
- **Neoepitope Presentation Analysis**  
- Support for RNA-seq data integration  

---

## **References**  📖💡
- GATK: [GATK Documentation](https://gatk.broadinstitute.org)  
- CNVkit: [CNVkit Docs](https://cnvkit.readthedocs.io)
- Sequenza:[Sequenza Documentation](https://bitbucket.org/sequenzatools/sequenza/src/chemins/)
- Control-FREEC: [FREEC](http://boevalab.com/FREEC/)  
- VEP: [VEP Documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)  
- SigProfiler: [SigProfilerExtractor](https://www.mathworks.com/help/bioinfo/ug/sigprofiler-extractor.html)
- Next-Flow-Nfcore: [Nfcore-piplines](https://nf-co.re/join)

---

## **Contact**  📞📧🤝💬
For inquiries or issues, contact me at: **chumawinnie@gmail.com**  

---  
