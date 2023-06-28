version = "1.1.0"
container_url = "nfontrodona/deseq2_wrapper:${version}"

params.formula = "~ condition"
params.comparison = "condition-TEST-CTRL"
params.gene_expression_threshold = 2
params.basemean_threshold = 0
params.lfc_threshold = 0
params.deseq2_out = ""
process deseq2_analysis {
  container = "${container_url}"
  label "small_mem_mono_cpus"

  if (params.deseq2_out != "") {
    publishDir "results/${params.deseq2_out}", mode: 'copy'
  }

  input:
    path design
    path filter
    path gene_name
    path count_file

  output:
    path "*.{html,pdf}", emit: figures
    path "Differential_expression_*_full.txt", emit: results
    path "*.csv", emit: counts
    path "*.txt", emit: statistics


  script:
    def pfilter = filter.name != 'NONE' ? "-f $filter" : ''
    def pgene = gene_name.name != 'EMPTY' ? "-c $gene_name" : ''
"""
Rscript /deseq2-wrapper.R ${design} "." "${params.formula}" "${params.comparison}" ${pfilter} ${pgene} -g ${params.gene_expression_threshold} -b ${params.basemean_threshold} -l ${params.lfc_threshold}
"""
}