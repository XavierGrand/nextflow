version = "1.0"
container_url = "xgrand/rtranger:${version}"

process porechop {
    container = "${container_url}"
    label "big_mem_mono_cpus"
    tag "$file_id"
    if (params.rtranger_out != "") {
    publishDir "results/${params.rtranger_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gtf)

  output:
    tuple val(file_id), path("*.bed"), emit: bed
  script:
"""
rtranger.py --gtf ${gtf} --output . --length ${rt_length}
"""
}