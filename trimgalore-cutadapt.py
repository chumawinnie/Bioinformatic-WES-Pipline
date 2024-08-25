import os
import yaml

# Read sample accession numbers from a file
with open("Access_id.txt") as file:
    accessions = file.read().split()  # Create a list of accessions

# Directory paths
inst_dir = os.path.expanduser("~/TrimGalore-0.6.7/trim_galore")
input_dir = "/home/obiorach/merged-whole-Exon-seq"
output_dir = "/home/obiorach/merged-whole-Exon-seq/trim_reads"
option_fastqc = "--fastqc"  # To use FastQC after trimming

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate over the list of accessions
for accession in accessions:
    # Define the paired input files
    input_file_1 = os.path.join(input_dir, f"{accession}_1.fastq.gz")
    input_file_2 = os.path.join(input_dir, f"{accession}_2.fastq.gz")

    # Determine sample type (Normal or Tumor) based on accession name
    sample_type = "Normal" if "Normal" in accession else "Tumor"

    if os.path.exists(input_file_1) and os.path.exists(input_file_2):
        # Construct the command for Trim Galore!
        command = (
            f"{inst_dir} "
            f"--paired "
            f"{input_file_1} {input_file_2} "
            f"-o {output_dir} "
            f"{option_fastqc} "
            f"-j 15"
        )

        print(f"Running {sample_type} sample: {command}")
        os.system(command)
        print(f"Finished running {sample_type} sample: {accession}")

        # Create a configuration dictionary
        config = {
            "accession": accession,
            "input_files": {
                "file_1": input_file_1,
                "file_2": input_file_2
            },
            "output_dir": output_dir,
            "options": {
                "fastqc": option_fastqc,
                "threads": 15
            },
            "sample_type": sample_type
        }
        
        # Write the configuration to a YAML file
        config_file_path = os.path.join(output_dir, f"{accession}_config.yml")
        with open(config_file_path, 'w') as config_file:
            yaml.dump(config, config_file, default_flow_style=False)

    else:
        print(f"Files {input_file_1} or {input_file_2} do not exist. Skipping this accession.")

print("All done")
print("This script is created by Chuma Winner Obiora, a Bioinformatician")
