# Antibody NGS Pipeline

Bulk antibody sequence preprocessing, annotation with abstar, import to MongoDB.  
This is based on the AbStar analysis pipeline: (www.github.com/briney/abstar)

### install  
`pip install antibody-ngs-pipeline`


### use  

To run antibody_ngs_pipeline:  
`antibody_ngs_pipeline`

To run antibody_ngs_pipeline with FASTQC report on raw data:  
`antibody_ngs_pipeline -f`
  
To run antibody_ngs_pipeline with adapter trimming by CutAdapt:  
`antibody_ngs_pipeline -t <path-to-adapters.fasta>`

To run antibody_ngs_pipeline with quality trimming by sickle:  
`antibody_ngs_pipeline -q`

To run antibody_ngs_pipeline with adapter trimming by CutAdapt, quality trimming 
with sickle and get a FASTQC report on both raw data and processed data:  
`antibody_ngs_pipeline -f -q -t <path-to-adapters.fasta>`




### requirements  
Python 3.5+  
abstar  
abutils  
cutadapt  
pymongo  
  

Downloading from BaseSpace requires Basemount: https://basemount.basespace.illumina.com/  
Quality trimming requires sickle: https://github.com/najoshi/sickle  
FastQC report requires FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
Merging paired sequences requires PANDAseq: https://github.com/neufeld/pandaseq  
batch_mongoimport (from Abstar) requires MongoDB: http://www.mongodb.org/  
