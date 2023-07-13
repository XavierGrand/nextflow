version = "1.0"
container_url = "xgrand/fusion_parser:${version}"

params.fusion_out = ""
process parsefusion{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "parse fusion"
  if (params.fusion_out != "") {
    publishDir "results/${params.fusion_out}", mode: 'copy'
  }

  input:
  path(fusions)
  path(count)
  path(design)

  output:
  path("*.tsv"), emit: dfg

  script:
  memory = "${task.memory}" - ~/\s*GB/

"""
Rscript fusion_parser.R --design ${design} --threads ${task.cpus} --memory ${memory}
"""
}