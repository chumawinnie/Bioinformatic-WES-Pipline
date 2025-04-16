import gzip
import os

# Define paths with expansion of '~'
vcf_file = os.path.expanduser("~/test-work-sarek/20results-folder/annotation/mutect2/Tumour_vs_Normal/Tumour_vs_Normal.mutect2.filtered_VEP.ann.vcf.gz")
bed_file = os.path.expanduser("~/whole-Exon-single-seq/ref-genome/exom_targets.bed/HyperExomeV2_capture_targets.hg19.bed")

# Define consequence terms considered relevant for TMB
tmb_effects = {
    "missense_variant", "stop_gained", "stop_lost",
    "start_lost", "frameshift_variant", "inframe_insertion",
    "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant"
}

# --- Step 1: Compute total coding region size from BED ---
def get_coding_region_size(bed_path):
    total_bases = 0
    with open(bed_path, 'r') as bed_file:
        for line in bed_file:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split('\t')
                start = int(parts[1])
                end = int(parts[2])
                total_bases += end - start
    return total_bases / 1e6  # Convert to megabases

# --- Step 2: Count non-synonymous mutations with FILTER=PASS ---
def count_tmb_mutations(vcf_path, effects_set):
    mutation_count = 0
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            filter_status = columns[6]
            if filter_status != "PASS":
                continue
            info = columns[7]
            csq_fields = [x for x in info.split(";") if x.startswith("CSQ=")]
            if csq_fields:
                csq_values = csq_fields[0].replace("CSQ=", "").split(",")
                for annotation in csq_values:
                    fields = annotation.split("|")
                    if len(fields) > 1 and fields[1] in effects_set:
                        mutation_count += 1
                        break  # Count once per variant
    return mutation_count

# --- Run pipeline ---
coding_region_mb = get_coding_region_size(bed_file)
mutation_count = count_tmb_mutations(vcf_file, tmb_effects)
tmb_value = mutation_count / coding_region_mb

# --- Output ---
print(f"🧬 Total non-synonymous mutations (FILTER=PASS): {mutation_count}") # subsetting the pass region from the vcf-file
print(f"📏 Coding region size: {coding_region_mb:.2f} Mb") # total coding region permegaBase from Bedfiles
print(f"✅ TMB-Score: {tmb_value:.2f} mutations/Mb")
print(f"🎉🎊Successfully done!.......\n🎖️ This Script was developed by chuma-winner, Bioinformatician @Uni-Augsburg Bioinformatic core-facility ")
