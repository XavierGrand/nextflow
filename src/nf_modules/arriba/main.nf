version = "2.3.0"
// container_url = "lbmc/arriba:${version}"
container_url = "xgrand/arriba:${version}"

params.arriba_options = "-f blacklist -E 1 -R 5000 -A 15 -M 2 -U 1000"
process arriba{
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"
  if (params.arriba_out != "") {
    publishDir "results/${params.arriba_out}", mode: 'copy'
  }

  input:
  tuple val(bam_id), path(bam), path(bai)
  path(gtf)
  tuple val(genome_id), path(genome)

  output:
  tuple val(bam_id), path("${bam_id}_fusions.tsv"), path("${bam_id}_fusions.discarded.tsv"), emit: fusions

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
  tuple val(fusion_id), path(fusions), path(discarded)
  tuple val(bam_id), path(bam), path(bai)
  path(gtf)
  path(cyto)

  output:
  path("${fusion_id}_fusions.pdf")

  script:
"""
Rscript /usr/local/bin/arriba_v${version}/draw_fusions.R \
    --fusions=${fusions} \
    --alignments=$bam \
    --output=${fusion_id}_fusions.pdf \
    --annotation=${gtf} \
    --cytobands=${cyto}
"""
}