version = "2.3.0"
// container_url = "lbmc/arriba:${version}"
container_url = "xgrand/arriba:${version}"

params.arriba_options = "-f blacklist"
process arriba{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${bam_id}"
  if (params.arriba_out != "") {
    publishDir "results/${params.arriba_out}", mode: 'copy'
  }

  input:
  tuple val(bam_id), path(bam)
  path(gtf)
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

process draw_fusions{
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"
  if (params.draw_fusions_out != "") {
    publishDir "results/${params.draw_fusions_out}", mode: 'copy'
  }

  input:
  tuple val(fusion_id), path(fusions)
  tuple val(bam_id), path(bam)
  path(gtf)

  output:
  tuple val(fusion_id), path("*.pdf"), emit: drawn_fusions

  script:
"""
./draw_fusions.R \
    --fusions=${fusions} \
    --alignments=${bam} \
    --output=${fusion_id}_fusions.pdf \
    --annotation=${gtf}
"""
}