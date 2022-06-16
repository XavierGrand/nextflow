version = "2.3.0"
// container_url = "lbmc/arriba:${version}"
container_url = "xgrand/arriba:${version}"

process arriba{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}"
  if (params.arriba_out != "") {
    publishDir "results/${params.arriba_out}", mode: 'copy'
  }

  input:
  tuple val(fasta_id), path(bam)
  tuple val(file_id), path(gtf)
  tuple val(file_id), path(genome)

  output:
  tuple val(file_id), path("fusions.tsv"), emit: fusions
  tuple val(file_id), path("fusions.discarded.tsv"), emit: discarded

  script:
"""
arriba -x ${bam} \
       -g ${gtf} \
       -a ${genome} \
       -o fusions.tsv \
       -O fusions.discarded.tsv
"""
}