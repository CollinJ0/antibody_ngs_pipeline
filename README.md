# Antibody NGS Pipeline

Bulk antibody sequence preprocessing, annotation with abstar, upload to MongoDB and S3

### install  
`pip install antibody-ngs-pipeline`


### use  

To run antibody_ngs_pipeline:  
`antibody_ngs_pipeline`

To run antibody_ngs_pipeline with FASTQC report on raw data:  
`antibody_ngs_pipeline -f`
  
To run antibody_ngs_pipeline with adapter trimming by CutAdapt:  
`antibody_ngs_pipeline -t <path-to-adapters.fasta>`

To run antibody_ngs_pipeline with adapter trimming by CutAdapt  
and FASTQC report on both raw data and adapter trimmed data:  
`antibody_ngs_pipeline -f -t <path-to-adapters.fasta>`
