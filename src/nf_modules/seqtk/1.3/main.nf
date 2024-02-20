version = "1.3"
container_url = "xgrand/seqtk:${version}"

params.sample_fastq_out = ""
params.reads_number = 100
process sample_fastq {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${barcode}"
  if (params.sample_fastq_out != "") {
    publishDir "results/${params.sample_fastq_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(fastq)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_sampled.fasta"), emit: sampled_fastq

  script:
    """
    mkdir ${barcode}
    seqtk sample -s100 ${fastq} ${params.reads_number} > ${barcode}/${barcode}_sampled.fastq
    seqtk seq -a ${barcode}/${barcode}_sampled.fastq > ${barcode}/${barcode}_sampled.fasta
    """
}
