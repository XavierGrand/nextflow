version = "0.6.0"
container_url = "lbmc/rasusa:${version}"

include { index_fasta } from "./../samtools/main.nf"

params.sample_fastq = ""
params.sample_fastq_coverage = ""
params.sample_fastq_size = ""
params.sample_fastq_out = ""
workflow sample_fastq {
  take:
  fastq
  fasta

  main:
  if (params.sample_fastq_coverage == "" && params.sample_fastq_size == ""){
    fastq
      .set{ final_fastq }
  } else {
    index_fasta(fasta)
    sub_sample_fastq(fastq, index_fasta.out.index)
    sub_sample_fastq.out.fastq
      .set{ final_fastq }
  }

  emit:
  fastq = final_fastq

}

process sub_sample_fastq {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.sample_fastq_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(index_id), path(idx)

  output:
    tuple val(file_id), path("sub_*.fastq.gz"), emit: fastq

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  sample_option = "-c " + params.sample_fastq_coverage
  if (params.sample_fastq_size != ""){
    sample_option = "-b " + params.sample_fastq_size
  }

  if (fastq.size() == 2)
"""
rasusa \
  -i ${fastq[0]} ${fastq[1]} \
  -g ${idx} \
  ${sample_option} \
  -o sub_${fastq[0].simpleName}.fastq.gz sub_${fastq[1].simpleName}.fastq.gz
"""
  else
"""
rasusa \
  -i ${fastq} \
  -g ${idx} \
  ${sample_option} \
  -o sub_${fastq.simpleName}.fastq.gz
"""
}