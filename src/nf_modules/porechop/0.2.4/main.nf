version = "0.2.4"
container_url = "xgrand/porechop:${version}"

params.porechop_out = ""
process porechop {
    container = "${container_url}"
    label "big_mem_multi_cpus"
    tag "$barcode"
    if (params.porechop_out != "") {
      publishDir "results/${params.porechop_out}", mode: 'copy'
    }

  input:
    tuple val(barcode), path(fastq)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_merged_porechoped.fastq.gz"), emit: porechoped_fastq
  script:
"""
mkdir ${barcode}
cd ${barcode}/
porechop --input ../${fastq} -o ${barcode}_merged_porechoped.fastq.gz --threads ${task.cpus}
"""
}