import os
import yaml

# Read sample accession numbers from a file
with open("Access_id.txt") as file:
    accessions = file.read().split()  # Create a list of accessions

# Directory paths
genome_index = "/home/obiorach/whole-Exon-single-seq/ref-genome/index/hg19.fa" #absolute part 4rm another folder
known_sites_snps = "/home/obiorach/whole-Exon-single-seq/ref-genome/known_sites.vcf/dbsnp_138.hg19.vcf" #absolute part 4rm another folder
known_sites_indels = "/home/obiorach/whole-Exon-single-seq/ref-genome/known_indels.vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" #absolute part 4rm another folder
exome_targets = "/home/obiorach/whole-Exon-single-seq/ref-genome/exom_targets.bed/HyperExomeV2_capture_targets.hg19.bed" #absolute part 4rm another folder
output_dir = "/home/obiorach/merged-whole-Exon-seq/mapping/"
output_bam = os.path.join(output_dir, "bam_sorted/")
output_vcf_dir = os.path.join(output_dir, "vcf/")
output_filtered_vcf_dir = os.path.join(output_dir, "filtered_vcf/")
output_recal_data_dir = os.path.join(output_dir, "recal_data/")
output_dedup_bam_dir = os.path.join(output_dir, "bam_dedup/")
qc_dir = os.path.join(output_dir, "qc/")

# Path to Picard JAR
picard_jar = "/home/obiorach/miniconda3/envs/gatk4/share/gatk4-4.5.0.0-0/picard.jar"

# Number of threads to use
num_threads = 20

# Ensure output directories exist
for directory in [output_vcf_dir, output_filtered_vcf_dir, output_recal_data_dir, output_dedup_bam_dir, qc_dir]:
    os.makedirs(directory, exist_ok=True)

# Initialize variables for tumor and normal samples
tumor_sample = None
normal_sample = None

# Iterate over the list of accessions
for accession in accessions:
    # Determine sample type (Normal or Tumor) based on accession name
    sample_type = "Normal" if "Normal" in accession else "Tumor"

    if sample_type == "Tumor":
        tumor_sample = accession
    elif sample_type == "Normal":
        normal_sample = accession

    if tumor_sample and normal_sample:
        # Define the output files for tumor and normal
        tumor_bam_file = os.path.join(output_bam, f"{tumor_sample}_sorted_rg.bam")
        normal_bam_file = os.path.join(output_bam, f"{normal_sample}_sorted_rg.bam")

        # Check if files exist before running the commands
        if not (os.path.exists(tumor_bam_file) and os.path.exists(normal_bam_file)):
            print(f"Error: BAM files for {tumor_sample} or {normal_sample} do not exist.")
            continue

        # Index known sites and indels VCF files
        if not os.path.exists(known_sites_snps + ".idx"):
            print(f"Indexing {known_sites_snps}")
            os.system(f"gatk IndexFeatureFile -I {known_sites_snps}")
        if not os.path.exists(known_sites_indels + ".idx"):
            print(f"Indexing {known_sites_indels}")
            os.system(f"gatk IndexFeatureFile -I {known_sites_indels}")

        # Remove duplicates for tumor
        tumor_dedup_bam_file = os.path.join(output_dedup_bam_dir, f"{tumor_sample}_dedup.bam")
        tumor_metrics_file = os.path.join(output_dedup_bam_dir, f"{tumor_sample}_dedup_metrics.txt")
        mark_duplicates_command = (
            f"java -Xmx8g -jar {picard_jar} MarkDuplicates I={tumor_bam_file} O={tumor_dedup_bam_file} "
            f"M={tumor_metrics_file} REMOVE_DUPLICATES=true"
        )
        print(f"Running: {mark_duplicates_command}")
        os.system(mark_duplicates_command)

        # Remove duplicates for normal
        normal_dedup_bam_file = os.path.join(output_dedup_bam_dir, f"{normal_sample}_dedup.bam")
        normal_metrics_file = os.path.join(output_dedup_bam_dir, f"{normal_sample}_dedup_metrics.txt")
        mark_duplicates_command = (
            f"java -Xmx8g -jar {picard_jar} MarkDuplicates I={normal_bam_file} O={normal_dedup_bam_file} "
            f"M={normal_metrics_file} REMOVE_DUPLICATES=true"
        )
        print(f"Running: {mark_duplicates_command}")
        os.system(mark_duplicates_command)

        # **Added lines for indexing deduplicated BAM files before BQSR**
        print(f"Indexing deduplicated BAM for tumor sample: {tumor_dedup_bam_file}")
        os.system(f"samtools index {tumor_dedup_bam_file}")
        print(f"Indexing deduplicated BAM for normal sample: {normal_dedup_bam_file}")
        os.system(f"samtools index {normal_dedup_bam_file}")

        # QC before BQSR
        qc_tumor_pre_realign = os.path.join(qc_dir, f"{tumor_sample}_pre_realign.flagstat")
        qc_normal_pre_realign = os.path.join(qc_dir, f"{normal_sample}_pre_realign.flagstat")
        print(f"Running: samtools flagstat {tumor_dedup_bam_file} > {qc_tumor_pre_realign}")
        os.system(f"samtools flagstat {tumor_dedup_bam_file} > {qc_tumor_pre_realign}")
        print(f"Running: samtools flagstat {normal_dedup_bam_file} > {qc_normal_pre_realign}")
        os.system(f"samtools flagstat {normal_dedup_bam_file} > {qc_normal_pre_realign}")

        # Base recalibration for tumor
        tumor_recal_data_file = os.path.join(output_recal_data_dir, f"{tumor_sample}_recal_data.table")
        tumor_recal_bam_file = os.path.join(output_dedup_bam_dir, f"{tumor_sample}_recal.bam")
        base_recalibrator_command = (
            f"gatk BaseRecalibrator -I {tumor_dedup_bam_file} -R {genome_index} "
            f"--known-sites {known_sites_snps} --known-sites {known_sites_indels} -L {exome_targets} -O {tumor_recal_data_file}"
        )
        print(f"Running: {base_recalibrator_command}")
        os.system(base_recalibrator_command)

        apply_bqsr_command = (
            f"gatk ApplyBQSR -I {tumor_dedup_bam_file} -R {genome_index} "
            f"--bqsr-recal-file {tumor_recal_data_file} -O {tumor_recal_bam_file}"
        )
        print(f"Running: {apply_bqsr_command}")
        os.system(apply_bqsr_command)

        # Base recalibration for normal
        normal_recal_data_file = os.path.join(output_recal_data_dir, f"{normal_sample}_recal_data.table")
        normal_recal_bam_file = os.path.join(output_dedup_bam_dir, f"{normal_sample}_recal.bam")
        base_recalibrator_command = (
            f"gatk BaseRecalibrator -I {normal_dedup_bam_file} -R {genome_index} "
            f"--known-sites {known_sites_snps} --known-sites {known_sites_indels} -L {exome_targets} -O {normal_recal_data_file}"
        )
        print(f"Running: {base_recalibrator_command}")
        os.system(base_recalibrator_command)

        apply_bqsr_command = (
            f"gatk ApplyBQSR -I {normal_dedup_bam_file} -R {genome_index} "
            f"--bqsr-recal-file {normal_recal_data_file} -O {normal_recal_bam_file}"
        )
        print(f"Running: {apply_bqsr_command}")
        os.system(apply_bqsr_command)

        # Re-add read groups after recalibration
        picard_command_tumor = (
            f"java -jar {picard_jar} AddOrReplaceReadGroups I={tumor_recal_bam_file} "
            f"O={tumor_recal_bam_file.replace('_recal.bam', '_recal_with_rg.bam')} "
            f"RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={tumor_sample}"
        )
        print(f"Running: {picard_command_tumor}")
        os.system(picard_command_tumor)

        picard_command_normal = (
            f"java -jar {picard_jar} AddOrReplaceReadGroups I={normal_recal_bam_file} "
            f"O={normal_recal_bam_file.replace('_recal.bam', '_recal_with_rg.bam')} "
            f"RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={normal_sample}"
        )
        print(f"Running: {picard_command_normal}")
        os.system(picard_command_normal)

        # Reindex the BAM files with corrected read groups
        reindexed_tumor_bam_file = tumor_recal_bam_file.replace('_recal.bam', '_recal_with_rg.bam')
        reindexed_normal_bam_file = normal_recal_bam_file.replace('_recal.bam', '_recal_with_rg.bam')
        os.system(f"samtools index {reindexed_tumor_bam_file}")
        os.system(f"samtools index {reindexed_normal_bam_file}")

        # Define the output files for somatic variant calling
        output_vcf_file = os.path.join(output_vcf_dir, f"{tumor_sample}_somatic.vcf")
        output_filtered_vcf_file = os.path.join(output_filtered_vcf_dir, f"{tumor_sample}_somatic_filtered.vcf")

        # Somatic variant calling using GATK Mutect2
        gatk_command = (
            f"gatk Mutect2 -R {genome_index} -I {reindexed_tumor_bam_file} -I {reindexed_normal_bam_file} "
            f"-tumor {tumor_sample} -normal {normal_sample} -L {exome_targets} -O {output_vcf_file} --native-pair-hmm-threads {num_threads}"
        )
        print(f"Running: {gatk_command}")
        os.system(gatk_command)

        # Filter somatic variants using FilterMutectCalls
        filter_command = (
            f"gatk FilterMutectCalls -V {output_vcf_file} -R {genome_index} -O {output_filtered_vcf_file}"
        )
        print(f"Running: {filter_command}")
        os.system(filter_command)

        # Create a configuration dictionary
        config = {
            "tumor_sample": tumor_sample,
            "normal_sample": normal_sample,
            "bam_files": {
                "tumor_bam": tumor_bam_file,
                "normal_bam": normal_bam_file,
                "tumor_dedup_bam": tumor_dedup_bam_file,
                "normal_dedup_bam": normal_dedup_bam_file,
                "tumor_recal_bam": tumor_recal_bam_file,
                "normal_recal_bam": normal_recal_bam_file,
                "reindexed_tumor_bam": reindexed_tumor_bam_file,
                "reindexed_normal_bam": reindexed_normal_bam_file
            },
            "vcf_files": {
                "output_vcf_file": output_vcf_file,
                "filtered_vcf_file": output_filtered_vcf_file
            },
            "recal_data_files": {
                "tumor_recal_data": tumor_recal_data_file,
                "normal_recal_data": normal_recal_data_file
            },
            "known_sites": {
                "snps": known_sites_snps,
                "indels": known_sites_indels
            },
            "exome_targets": exome_targets,
            "num_threads": num_threads
        }

        # Write the configuration to a YAML file
        config_file_path = os.path.join(output_dir, f"{tumor_sample}_config.yml")
        with open(config_file_path, 'w') as config_file:
            yaml.dump(config, config_file, default_flow_style=False)

        # Reset for next pair
        tumor_sample = None
        normal_sample = None

print("Somatic Exon variant calling and filtering completed.")
print("Created by Chuma Winner Obiora, a bioinformatician")
