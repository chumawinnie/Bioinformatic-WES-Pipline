nextflow run nf-core/hlatyping -r 2.0.0 -profile docker \
  --input /home/obiorach/test-work-sarek/3results-folder/sample3-HLA-analysis.csv \
  --outdir /home/obiorach/test-work-sarek/3results-folder/HLA3-analysis-results \
  --genome hg19 \
  --max_cpus 20 \
  --max_memory '30 GB' \
  -c custom.HLA-analysis.config

###here is the content of the config file 
process {
    withName: ".*" {
        time = 72.h                                                                               
    }
}
