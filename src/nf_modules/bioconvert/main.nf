version = "0.4.0"
container_url = "lbmc/bioconvert:${version}"
params.bigwig_to_wig = ""
params.bigwig_to_wig_out = ""
process bigwig_to_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id) path(bw)

  output:
  tuple val(file_id), path("*.wig"), emit: wig

  script:
"""
bioconvert bigwig2wiggle ${bw} ${bw.simpleName}.wig
"""
}

params.bigwig2_to_wig2 = ""
params.bigwig2_to_wig2_out = ""
process bigwig2_to_wig2 {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bw_a), path(bw_b)

  output:
  tuple val(file_id), path("${bw_a.simpleName}.wig"), path("${bw_b.simpleName}.wig"), emit: wig

  script:
"""
bioconvert bigwig2wiggle ${bw_a} ${bw_a.simpleName}.wig
bioconvert bigwig2wiggle ${bw_b} ${bw_b.simpleName}.wig
"""
}