What you will find in this Repo : ready to use code snippets for processing RNAseq ( bulk )  datasets using Nextflow workflow. 


  #### Git folder files and description 
      Nextflow-Bulk-RNAseq-workflow
      ├─ README.md - READ ME file
      ├─ main.nf   - nextflow workflow 
      ├─ nextflow.config - execution enviormemt 
      ├─ conf/ - global parameters
      │  ├─ base.config 
      │  ├─ slurm.config
      │  ├─ aws.config
      │  └─ gcp.config
      ├─ assets/ - example sample format and multiqc parameters 
      │  ├─ samplesheet.csv
      │  └─ multiqc_config.yaml
      ├─ .gitignore - list of files / folders not to be tracked

#### RNA-seq Nextflow (DSL2)
    bulk RNA-seq primary analysis 
                            |_adapter trimming 
                              |_fastQC ( pre and post adapter trim ) 
                                |_Alignment ( STAR ) 
                                  |_quantification (featurecount) 
                                    |_multiQC
        

#### Requirements 
    Nextflow >= 23.10
    Containers: Docker or Singularity/Apptainer
    References: genome FASTA and matched GTF
    Cloud: AWS/GCP credentials configured if using those profiles

#### nextflow command
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




