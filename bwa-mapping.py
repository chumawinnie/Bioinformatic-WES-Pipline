import os
import yaml

# Read sample accession numbers from a file
with open("Access_id.txt") as file:
    accessions = file.read().split()  # Create a list of accessions

# Directory paths
inst_dir = "bwa"  # Assuming BWA is installed and accessible from PATH
genome_index = "/home/obiorach/whole-Exon-single-seq/ref-genome/index/hg19.fa"  # Path to BWA human-genome GRCh38 index which is in another folder!(absolute path)
input_dir = "/home/obiorach/merged-whole-Exon-seq/trim_reads/"
output_dir = "/home/obiorach/merged-whole-Exon-seq/mapping/"
output_sam_dir = os.path.join(output_dir, "sam/")
output_sorted_bam_dir = os.path.join(output_dir, "bam_sorted/")
output_bam_dir = os.path.join(output_dir, "bam/")
output_unaligned = os.path.join(output_dir, "unaligned/")
reports_dir = os.path.join(output_dir, "reports/")
cores = 20  # Number of threads for BWA

# Specify the path to the GATK executable
gatk_path = "/home/obiorach/miniconda3/envs/gatk4/bin/gatk"

# Ensure output directories exist
for directory in [output_dir, output_sam_dir, output_sorted_bam_dir, output_bam_dir, output_unaligned, reports_dir]:
    os.makedirs(directory, exist_ok=True)

# Iterate over the list of accessions
for accession in accessions:
    # Define the paired input files
    input_file_1 = os.path.join(input_dir, f"{accession}_1_val_1.fq.gz")
    input_file_2 = os.path.join(input_dir, f"{accession}_2_val_2.fq.gz")

    # Determine sample type (Normal or Tumor) based on accession name
    sample_type = "Normal" if "Normal" in accession else "Tumor"

    # Define the output files
    output_sam_file = os.path.join(output_sam_dir, f"{accession}.sam")
    output_sorted_bam_file = os.path.join(output_sorted_bam_dir, f"{accession}_sorted.bam")
    output_bam_file = os.path.join(output_bam_dir, f"{accession}.bam")
    output_rg_bam_file = os.path.join(output_sorted_bam_dir, f"{accession}_sorted_rg.bam")

    if os.path.exists(input_file_1) and os.path.exists(input_file_2):
        # Map reads using BWA
        bwa_command = (
            f"{inst_dir} mem -t {cores} {genome_index} {input_file_1} {input_file_2} > {output_sam_file}"
        )
        print(f"Running {sample_type} sample BWA mapping: {bwa_command}")
        os.system(bwa_command)

        # Sort the SAM file and convert to BAM
        sort_command = f"samtools sort -@ {cores} -o {output_sorted_bam_file} {output_sam_file}"
        print(f"Running {sample_type} sample SAM sorting: {sort_command}")
        os.system(sort_command)

        # Add Read Groups using GATK
        rgid = "group1"
        rglb = "lib1"
        rgpl = "illumina"
        rgpu = "unit1"
        rgsample = accession  # You can customize this based on your needs

        add_rg_command = (
            f"{gatk_path} AddOrReplaceReadGroups "
            f"I={output_sorted_bam_file} "
            f"O={output_rg_bam_file} "
            f"RGID={rgid} "
            f"RGLB={rglb} "
            f"RGPL={rgpl} "
            f"RGPU={rgpu} "
            f"RGSM={rgsample}"
        )
        print(f"Adding Read Groups to {sample_type} sample: {add_rg_command}")
        os.system(add_rg_command)

        # Index the BAM file with read groups
        index_rg_command = f"samtools index {output_rg_bam_file}"
        print(f"Indexing BAM file with Read Groups for {sample_type} sample: {index_rg_command}")
        os.system(index_rg_command)

        # Create a configuration dictionary
        config = {
            "accession": accession,
            "input_files": {
                "file_1": input_file_1,
                "file_2": input_file_2
            },
            "output_files": {
                "sam_file": output_sam_file,
                "sorted_bam_file": output_sorted_bam_file,
                "rg_bam_file": output_rg_bam_file
            },
            "output_dir": output_dir,
            "options": {
                "bwa_threads": cores,
                "fastqc": "N/A",  # No FastQC in this script
                "gatk_path": gatk_path,
                "rg_info": {
                    "RGID": rgid,
                    "RGLB": rglb,
                    "RGPL": rgpl,
                    "RGPU": rgpu,
                    "RGSM": rgsample
                }
            },
            "sample_type": sample_type
        }

        # Write the configuration to a YAML file
        config_file_path = os.path.join(output_dir, f"{accession}_config.yml")
        with open(config_file_path, 'w') as config_file:
            yaml.dump(config, config_file, default_flow_style=False)

    else:
        print(f"Files {input_file_1} or {input_file_2} do not exist. Skipping this accession.")

print("Mapping, sorting, adding read groups, and indexing completed.")
print("All done")
print("This script was created by Chuma Winner Obiora, a Bioinformatician @ Uni-Augsburg")  
print("Supervised and Approved by:")
print("Dr. Jan Meier-Kolthoff: Head of Bioinformatics Core Facility @ Uni-Augsburg")
print("Prof. Dr. Matthias Schlesner: Chair of Biomedical Informatics, Data Mining and Data Analytics @ University of Augsburg")
