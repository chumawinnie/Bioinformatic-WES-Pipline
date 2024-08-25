import os
import glob

# All files ending with .gz/fastq
files = glob.glob("/home/obiorach/merged-whole-Exon-seq/*.fastq.gz")

inst_dir = "/usr/bin/fastqc"
input_dir = "/home/obiorach/merged-whole-Exon-seq"
output_dir = "/home/obiorach/merged-whole-Exon-seq/fastqc-output"

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Loop through the files to check the quality of the files
for x in files:
    command = f"{inst_dir} {x} -o {output_dir}"
    print(f"{x} ........ Running the first quality control of the sequencing data")
    os.system(command)

print("All done")
print("This script was created by Chuma Winner Obiora, a Bioinformatician @ Uni-Augsburg")
print("Supervised and Approved by:")
print("Dr. Jan Meier-Kolthoff: Head of Bioinformatics Core Facility @ Uni-Augsburg")
print("Prof. Dr. Matthias Schlesner: Chair of Biomedical Informatics, Data Mining and Data Analytics @ University of Augsburg")
