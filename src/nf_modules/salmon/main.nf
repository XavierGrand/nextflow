version = "1.8.0"
container_url = "lbmc/salmon:${version}"

process quantify {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.salmon_out != "") {
    publishDir "results/${params.salmon_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.sf"), emit: quant
  script:
"""
salmon quant -l A --noErrorModel -t XXXXXXXXXX -a ${bam} -p 4 -o ${params.salmon_out}
"""
}