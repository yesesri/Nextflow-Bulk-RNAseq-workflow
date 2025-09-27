# RNA-seq Nextflow (DSL2)
bulk RNA-seq primary analysis 
                            |_adapter trimming 
                              |_fastQC
                                |_Alignment ( STAR ) 
                                  |_quantification ( Salmon ) 
                                    |_multiQC
        

# Requirements 
    1) Install Nextflow (>=23.10) and Docker/Singularity
     2) Clone repo and cd into it

# nextflow command
    nextflow run . \
    --samplesheet assets/samplesheet.csv \
    --fasta /path/to/genome.fa \
    --gtf   /path/to/genes.gtf \
    -profile docker,local
