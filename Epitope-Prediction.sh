nextflow run nf-core/epitopeprediction -r 2.3.1 -profile docker \
  --input /home/obiorach/test-work-sarek/3results-folder/sample3-Epitope-prediction.csv \
  --outdir /home/obiorach/test-work-sarek/3results-folder/Epitope-prediction-results \
  --tools netmhcpan-4.1 \
  --netmhcpan_path /home/obiorach/netMHCpan-4.1 \
  --max_cpus 20 \
  --max_memory '30 GB' \
  -c /home/obiorach/test-work-sarek/3results-folder/custom.Epitope-prediction.config
