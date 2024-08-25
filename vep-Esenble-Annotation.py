import os  # Import the os module to interact with the operating system

# Directory paths
output_filtered_vcf_dir = "/home/obiorach/merged-whole-Exon-seq/mapping/filtered_vcf/"  # Path to the directory containing filtered VCF files
vep_output_dir = "/home/obiorach/merged-whole-Exon-seq/mapping/vep_annotations/"  # Path to the directory where annotated VCF files will be saved
vep_cache_dir = "/home/obiorach/vep_cache"  # Path to the directory where the VEP cache is stored
vep_executable = "/home/obiorach/ensembl-vep/vep"  # Path to the VEP script

# Ensure the VEP output directory exists
os.makedirs(vep_output_dir, exist_ok=True)  # Create the VEP output directory if it doesn't exist

# Number of threads to use
num_threads = 20  # Update this to the number of CPU cores you want to use

# Iterate over all VCF files in the filtered VCF directory
for filename in os.listdir(output_filtered_vcf_dir):  # List all files in the filtered VCF directory
    if filename.endswith(".vcf"):  # Check if the file has a .vcf extension
        input_vcf_path = os.path.join(output_filtered_vcf_dir, filename)  # Construct the full path to the input VCF file
        output_vcf_path = os.path.join(vep_output_dir, filename.replace(".vcf", "_annotated.vcf"))  # Construct the full path to the output annotated VCF file

        # Construct the VEP command
        vep_command = (
            f"{vep_executable} -i {input_vcf_path} -o {output_vcf_path} --vcf --cache --dir_cache {vep_cache_dir} "  # Basic command structure with input and output paths
            f"--assembly GRCh37 --symbol --canonical --protein --biotype --uniprot --hgvs --fork {num_threads}"  # Additional options for VEP to include various annotation details
        )

        # Print and run the VEP command
        print(f"Running: {vep_command}")  # Print the VEP command to the console for verification
        os.system(vep_command)  # Execute the VEP command using the os.system() function

print("VEP annotation completed.")  # Print a message indicating the annotation process is complete
print("This script was created by Chuma Winner Obiora, a Bioinformatician @ Uni-Augsburg")
print("Supervised and Approved by:")
print("Dr. Jan Meier-Kolthoff: Head of Bioinformatics Core Facility @ Uni-Augsburg")
print("Prof. Dr. Matthias Schlesner: Chair of Biomedical Informatics, Data Mining and Data Analytics @ University of Augsburg")

