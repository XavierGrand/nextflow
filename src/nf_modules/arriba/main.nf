version = "2.3.0"
// container_url = "lbmc/arriba:${version}"
container_url = "xgrand/arriba:${version}"

params.arriba_options = "-f blacklist"
process arriba{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}"
  if (params.arriba_out != "") {
    publishDir "results/${params.arriba_out}", mode: 'copy'
  }

  input:
  tuple val(bam_id), path(bam)
  tuple val(gtf_id), path(gtf)
  tuple val(genome_id), path(genome)

  output:
  tuple val(bam_id), path("*fusions.tsv"), emit: fusions
  tuple val(bam_id), path("*fusions.discarded.tsv"), emit: discarded

  script:
"""
arriba -x ${bam} \
       -g ${gtf} \
       -a ${genome} \
       -o ${bam_id}_fusions.tsv \
       -O ${bam_id}_fusions.discarded.tsv \
       ${params.arriba_options}
"""
}