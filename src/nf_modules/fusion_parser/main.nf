version = "1.0"
container_url = "xgrand/fusion_parser:${version}"

process parsefusion{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${bam_id}"
  if (params.fusion_out != "") {
    publishDir "results/${params.fusion_out}", mode: 'copy'
  }

  input:
  path(fusion)
  path(count)
  path(design)

  output:
  tuple val(bam_id), path("*.tsv"), emit: dfg

  script:
  memory = "${task.memory}" - ~/\s*GB/
"""
Rscript fusion_parser.R --fusion ${fusion} --count ${count} --design ${design} --threads ${task.cpus} --memory ${memory}
"""
}