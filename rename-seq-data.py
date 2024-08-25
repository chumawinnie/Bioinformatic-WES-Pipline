import os

# Define the old and new filenames
old_filenames = [
    "3-N_R1_001.fastq.gz", 
    "3-N_R2_001.fastq.gz", 
    "3_R1_001.fastq.gz", 
    "3_R2_001.fastq.gz"
]

new_filenames = {
    "3_R1_001.fastq.gz": "3-R1_001-Tumour_1.fastq.gz",
    "3_R2_001.fastq.gz": "3-R2_001-Tumour_2.fastq.gz",
    "3-N_R1_001.fastq.gz": "3-N_R1_001-Normal_1.fastq.gz",
    "3-N_R2_001.fastq.gz": "3-N_R2_001-Normal_2.fastq.gz"
}

# Rename the files
for old_name in old_filenames:
    if old_name in new_filenames:
        new_name = new_filenames[old_name]
        os.rename(old_name, new_name)
        print(f"Renamed {old_name} to {new_name}")
    else:
        print(f"No new name defined for {old_name}")

