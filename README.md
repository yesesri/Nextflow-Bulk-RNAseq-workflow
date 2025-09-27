# RNA-seq Nextflow (DSL2)
bulk RNA-seq primary analysis 
                            |_adapter trimming 
                              |_fastQC
                                |_Alignment ( STAR ) 
                                  |_quantification ( Salmon ) 
                                    |_multiQC
        

# Requirements 
    Nextflow >= 23.10
    Containers: Docker or Singularity/Apptainer
    References: genome FASTA and matched GTF
    Cloud: AWS/GCP credentials configured if using those profiles

# nextflow command
    nextflow run . \
    --samplesheet assets/samplesheet.csv \
    --fasta /path/to/genome.fa \
    --gtf   /path/to/genes.gtf \
    -profile 

    note : 
      - assets/samplesheet.csv :  example samplesheet, follow same format. 
      -  profile options  
        # SLURM (HPC)
        nextflow run . \
          --samplesheet assets/samplesheet.csv \
           --fasta /proj/ref/genome.fa --gtf /proj/ref/genes.gtf \


         # AWS Batch (adjust aws.config & S3 paths)
          nextflow run . \
            --samplesheet s3://my-bucket/samplesheet.csv \
            --fasta s3://ref-bucket/genome.fa --gtf s3://ref-bucket/genes.gtf \
            --outdir s3://my-bucket/results \
            -profile aws

        # GCP Life Sciences (adjust gcp.config & GCS paths)
        nextflow run . \
          --samplesheet gs://my-bucket/samplesheet.csv \
          --fasta gs://ref-bucket/genome.fa --gtf gs://ref-bucket/genes.gtf \
          --outdir gs://my-bucket/results \
          -profile gcp

