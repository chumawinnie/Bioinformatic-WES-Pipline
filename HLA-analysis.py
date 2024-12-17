import os  # Import the os module to interact with the operating system

# Define the paths to the paired-end FASTQ files for tumor and normal samples
normal_r1 = "/home/obiorach/datasets/WES_testdata_from_DNPM_trial/sample_no_3/fastq_data_from_single_seq_run/3-N_R1_001.fastq.gz"
normal_r2 = "/home/obiorach/datasets/WES_testdata_from_DNPM_trial/sample_no_3/fastq_data_from_single_seq_run/3-N_R2_001.fastq.gz"
tumor_r1 = "/home/obiorach/datasets/WES_testdata_from_DNPM_trial/sample_no_3/fastq_data_from_single_seq_run/3_R1_001.fastq.gz"
tumor_r2 = "/home/obiorach/datasets/WES_testdata_from_DNPM_trial/sample_no_3/fastq_data_from_single_seq_run/3_R2_001.fastq.gz"

# Output directories for OptiType results
normal_output_dir = "/home/obiorach/test-work-sarek/optitype_output/normal-results"
tumor_output_dir = "/home/obiorach/test-work-sarek/optitype_output/Tumour-results"

# Ensure output directories exist
os.makedirs(normal_output_dir, exist_ok=True)
os.makedirs(tumor_output_dir, exist_ok=True)

# Command template for running OptiType
optitype_executable = "/home/obiorach/OptiType/OptiTypePipeline.py"  # Full path to OptiTypePipeline.py

# Run OptiType for the normal sample
normal_cmd = (
    f"python3 {optitype_executable} -i {normal_r1} {normal_r2} -r -o {normal_output_dir}"
)
print(f"Running OptiType for Normal sample: {normal_cmd}")
os.system(normal_cmd)  # Execute the command to run OptiType for normal sample

# Run OptiType for the tumor sample
tumor_cmd = (
    f"python3 {optitype_executable} -i {tumor_r1} {tumor_r2} -r -o {tumor_output_dir}"
)
print(f"Running OptiType for Tumor sample: {tumor_cmd}")
os.system(tumor_cmd)  # Execute the command to run OptiType for tumor sample

# Print a completion message
print("OptiType HLA typing completed for both Normal and Tumor samples.")
print("This script was created by Chuma Winner Obiora, a Bioinformatician @ Uni-Augsburg")
print("Supervised and Approved by:")
print("Dr. Jan Meier-Kolthoff: Head of Bioinformatics Core Facility @ Uni-Augsburg")
print("Prof. Dr. Matthias Schlesner: Chair of Biomedical Informatics, Data Mining and Data Analytics @ University of Augsburg")
