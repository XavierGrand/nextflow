version = "1.0.0"
container_url = "nfontrodona/topgo_wrapper:${version}"

params.top = 20
params.id = "symbol"
params.basemean_threshold = 0
params.lfc_threshold = 0
params.topgo_out = ""
process topgo_analysis {
  container = "${container_url}"
  label "big_mem_mono_cpus"

  if (params.topgo_out != "") {
    publishDir "results/${params.topgo_out}", mode: 'copy'
  }

  input:
    path results

  output:
    path "*.txt", emit: tables
    path "*.pdf", emit: figures


  script:
"""
Rscript /topGO-wrapper.R -d ${results} -i ${params.id} -b ${params.basemean_threshold} -l ${params.lfc_threshold}
"""
}