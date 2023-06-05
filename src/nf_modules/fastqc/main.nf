version = "0.11.5"
container_url = "lbmc/fastqc:${version}"

params.fastqc_fastq = ""
params.fastqc_fastq_out = ""
process fastqc_fastq {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.fastqc_fastq_out != "") {
    publishDir "results/${params.fastqc_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.{zip,html}"), emit: report

  script:
  if (reads.size() == 2)
  """
  fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
    ${params.fastqc_fastq} \
    ${reads[0]} ${reads[1]}
  """
  else
  """
    fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ ${params.fastqc_fastq} ${reads[0]}
  """
}

process fastqc_fastq_multi {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.fastqc_fastq_out != "") {
    publishDir "results/${params.fastqc_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.{zip,html}"), emit: report

  script:
  if (reads.size() == 2)
  """
  fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
    ${params.fastqc_fastq} \
    ${reads[0]} ${reads[1]}
  """
  else
  """
    fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ ${params.fastqc_fastq} ${reads[0]}
  """
}