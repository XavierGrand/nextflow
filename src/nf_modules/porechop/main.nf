version = "0.2.4"
container_url = "lbmc/porechop:${version}"

process porechop {
    container = "${container_url}"
    label "big_mem_multi_cpus"
    tag "$file_id"
    if (params.porechop_out != "") {
    publishDir "results/${params.porechop_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fatsq)

  output:
    tuple val(file_id), path("*_porechoped.fastq"), emit: porechoped_fastq
  script:
"""
porechop -i ${fastq} -o ${file_id}_porechoped.fastq --threads 4
"""
}