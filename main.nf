/* enables DSL2 version of Nextflow */
nextflow.enable.dsl=2

/* parameters |  user Inputs - fasta , gtf files and path to output dir */
params.samplesheet   = "assets/samplesheet.csv"
params.outdir        = "results"
params.fasta         = null
params.gtf           = null
params.star_index    = null
params.max_cpus      = 8
params.max_memory    = "32 GB"
params.max_time      = "24h"

/* Manifest - details about the workflow,will be displayed in log and report files  */
manifest {
  name            = 'bulk-rnaseq-nf-dsl2'
  author          = 'Your Name'
  description     = 'Bulk RNA-seq: trim, QC, STAR align, featureCounts, MultiQC'
  nextflowVersion = '>=23.10.0'
}

/* Inputs & Channels 
Basic sanity checks (fail early with clear messages) */

assert params.fasta : "Missing --fasta path to reference genome FASTA"
assert params.gtf   : "Missing --gtf path to annotation GTF"

/* Samples: expect header sample,fastq1,fastq2 */
samples_ch = Channel
  .fromPath(params.samplesheet)
  .splitCsv(header:true)
  .map { row ->
    assert row.sample && row.fastq1 && row.fastq2 : "CSV must include columns: sample,fastq1,fastq2"
    tuple( row.sample as String, [ file(row.fastq1), file(row.fastq2) ] )
  }

fasta_ch = Channel.value( file(params.fasta) )
gtf_ch   = Channel.value( file(params.gtf) )

/* STAR index: use provided path or build from FASTA+GTF */
star_index_ch = params.star_index ?
  Channel.value( file(params.star_index) ) :
  MAKE_STAR_INDEX(fasta_ch, gtf_ch).out.index_dir

/* Processes */

/* run fastqc on I/P files */
process FASTQC_RAW {
  tag { sample_id }
  publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
  cpus 2
  memory '2 GB'
  time '2h'
  container 'quay.io/biocontainers/fastqc:0.12.1--0'

  input:
    tuple val(sample_id), path(reads)

  output:
    path "*.zip",  emit: zips
    path "*.html", emit: html

  script:
  """
  fastqc -t ${task.cpus} -o ./ ${reads[0]} ${reads[1]}
  """
}

/* trim adapters */
process TRIM_FASTP {
  tag { sample_id }
  publishDir "${params.outdir}/fastp", mode: 'copy'
  cpus 4
  memory '8 GB'
  time '6h'
  container 'quay.io/biocontainers/fastp:0.23.4--h2e03b76_1'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_pairs
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

  script:
  """
  fastp \
    -i ${reads[0]} -I ${reads[1]} \
    -o ${sample_id}_R1.trim.fastq.gz -O ${sample_id}_R2.trim.fastq.gz \
    --detect_adapter_for_pe -w ${task.cpus} \
    -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json
  """
}

/* quality control check point - post adapter trimming */

process FASTQC_TRIMMED {
  tag { sample_id }
  publishDir "${params.outdir}/fastqc_trim", mode: 'copy'
  cpus 2
  memory '2 GB'
  time '2h'
  container 'quay.io/biocontainers/fastqc:0.12.1--0'

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    path "*.zip",  emit: zips
    path "*.html", emit: html

  script:
  """
  fastqc -t ${task.cpus} -o ./ ${r1} ${r2}
  """
}

/* create index file for alignment */

process MAKE_STAR_INDEX {
  tag "star_index"
  publishDir "${params.outdir}/star_index", mode: 'copy'
  cpus 8
  memory '30 GB'
  time '24h'
  container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_3'

  input:
    path fasta
    path gtf

  output:
    path "star_index", emit: index_dir

  script:
  """
  mkdir -p star_index
  STAR --runMode genomeGenerate \
       --runThreadN ${task.cpus} \
       --genomeDir star_index \
       --genomeFastaFiles ${fasta} \
       --sjdbGTFfile ${gtf} \
       --sjdbOverhang 100
  """
}

/* Alignment */ 

process STAR_ALIGN {
  tag { sample_id }
  publishDir "${params.outdir}/star", mode: 'copy'
  cpus 8
  memory '24 GB'
  time '24h'
  container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_3'

  input:
    tuple val(sample_id), path(r1), path(r2)
    path star_index

  output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}_Log.final.out", emit: log_final
    path "${sample_id}_Log.out",       emit: log_out

  script:
  """
  STAR --genomeDir ${star_index} \
       --readFilesIn ${r1} ${r2} \
       --readFilesCommand zcat \
       --runThreadN ${task.cpus} \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ${sample_id}_

  mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
  """
}

/* index alinged bam files */

process INDEX_BAM {
  tag { sample_id }
  publishDir "${params.outdir}/star", mode: 'copy'
  cpus 1
  memory '2 GB'
  time '2h'
  container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: bam_indexed

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  """
}

/* Qualtification form aligned BAMS ; can choose to use salmon on raw fastq's based on project requirement */
process FEATURECOUNTS {
  tag "featureCounts"
  publishDir "${params.outdir}/featurecounts", mode: 'copy'
  cpus 8
  memory '8 GB'
  time '8h'
  container 'quay.io/biocontainers/subread:2.0.6--he4a0461_1'

  input:
    path bams   // list of BAMs (collected)
    path gtf

  output:
    path "gene_counts.txt",         emit: counts
    path "gene_counts.txt.summary", emit: summary

  script:
  """
  featureCounts -T ${task.cpus} \
    -a ${gtf} -t exon -g gene_id \
    -o gene_counts.txt \
    ${bams.join(' ')}
  """
}

/* generate QC report on all smaples in project or study */

process MULTIQC {
  tag "multiqc"
  publishDir "${params.outdir}/multiqc", mode: 'copy'
  cpus 2
  memory '2 GB'
  time '2h'
  container 'quay.io/biocontainers/multiqc:1.20--pyhdfd78af_0'

  input:
    path(qc_inputs)

  output:
    path "multiqc_report.html"
  script:
  """
  multiqc -c ${projectDir}/assets/multiqc_config.yaml -o ./ .
  """
}

/* WORKFLOW -  call each process in order of processing */
workflow {
  fastqc_raw = FASTQC_RAW(samples_ch)
  trimmed = TRIM_FASTP(samples_ch)
  fastqc_trim = FASTQC_TRIMMED(trimmed.trimmed_pairs)
  aligned = STAR_ALIGN(trimmed.trimmed_pairs, star_index_ch)
  indexed = INDEX_BAM(aligned.bam)
  all_bams = indexed.bam_indexed.map { tup -> tup[1] }.collect()  // extract BAM paths
  featc = FEATURECOUNTS(all_bams, gtf_ch)

  qc_mix = Channel
    .merge(
      fastqc_raw.zips, fastqc_raw.html,
      trimmed.html, trimmed.json,
      fastqc_trim.zips, fastqc_trim.html,
      aligned.log_final, aligned.log_out,
      featc.counts, featc.summary
    ).collect()

  MULTIQC(qc_mix)
}




