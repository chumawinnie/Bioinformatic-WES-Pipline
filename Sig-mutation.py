import os
import shutil
from SigProfilerExtractor import sigpro as sig
from SigProfilerMatrixGenerator import install as genInstall

def run_signature_extraction():
    # Directory paths
    vcf_dir = "/home/obiorach/whole-Exon-single-seq/mapping/vcf"
    output_dir = "/home/obiorach/whole-Exon-single-seq/mapping/vcf/output-result"
    specific_vcf_file = os.path.join(vcf_dir, "3_R1_001-Tumour_somatic.vcf") 

    # Create a temporary directory to process only the specific VCF file
    temp_vcf_dir = "/home/obiorach/whole-Exon-single-seq/mapping/vcf/temp_vcf_dir"
    os.makedirs(temp_vcf_dir, exist_ok=True)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Check if the specific VCF file exists
    if os.path.exists(specific_vcf_file):
        # Copy the specific VCF file into the temporary directory
        shutil.copy(specific_vcf_file, temp_vcf_dir)

        # Install the GRCh37 genome if not already installed
        print("Installing GRCh37 genome if not already installed...")
        genInstall.install('GRCh37')

        # Parameters for signature extraction
        minimum_signatures = 1
        maximum_signatures = 5
        nmf_replicates = 100

        # Limit to 1 CPU core to avoid multiprocessing errors
        print(f"Running mutational signature extraction on {specific_vcf_file}")
        sig.sigProfilerExtractor(
            "vcf", 
            output_dir, 
            temp_vcf_dir,  # Use the temporary directory containing only the specific VCF file
            minimum_signatures=minimum_signatures, 
            maximum_signatures=maximum_signatures, 
            nmf_replicates=nmf_replicates,
            cpu=1  # Limit to 1 CPU core
        )

        print(f"Finished mutational signature extraction. Results saved to {output_dir}")

        # Optionally, clean up the temporary directory after the process
        shutil.rmtree(temp_vcf_dir)
    else:
        print(f"VCF file {specific_vcf_file} does not exist. Please check the path and try again.")

if __name__ == "__main__":
    run_signature_extraction()
    print("This script is created by Chuma Winner Obiora, a Bioinformatician")
