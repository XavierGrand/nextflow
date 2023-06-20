version = "1.0"
container_url = "xgrand/rtranger:${version}"

params.rt_length=10000
process rtranger {
    container = "${container_url}"
    label "big_mem_mono_cpus"
    tag "rt range computation"
    if (params.rtranger_out != "") {
    publishDir "results/${params.rtranger_out}", mode: 'copy'
  }

  input:
    path(gtf)

  output:
    path("readthrough_range.tsv"), emit: rtrange_tsv
    path("first_last_exons.tsv"), emit: fl_exon_tsv
    path("*.bed"), emit: bed
  script:
"""
python ./rtranger.py --gtf ${gtf} --output . --length ${params.rt_length}
"""
}