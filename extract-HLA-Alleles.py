import pandas as pd

# Load the OptiType result file
hla_result_path = '~/test-work-sarek/optitype_output/Tumour-results/2024_10_14_13_57_41/2024_10_14_13_57_41_result.tsv'
hla_df = pd.read_csv(hla_result_path, sep='\t')

# Extract the HLA alleles
hla_alleles = hla_df[['A1', 'A2', 'B1', 'B2', 'C1', 'C2']]
print(hla_alleles)

