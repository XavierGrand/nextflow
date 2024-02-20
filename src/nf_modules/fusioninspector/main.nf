version = "2.8.0"
// container_url = "lbmc/star-fusion:${version}"
container_url = "xgrand/fusioninspector:${version}"

params.fusioninspector_out = "inspected_fusion"
process star_fusion {

    container = "${container_url}"
    label "big_mem_multi_cpus"
    tag "$file_id"
    if (params.fusioninspector_out != "") {
        publishDir "results/${params.fusioninspector_out}", mode: 'copy'
    }

    input:
        path(genome_lib)
        tuple val(file_id), path(fastq)

    output:
        tuple val(file_id), path ("*")

    script:

"""
FusionInspector --fusions fusions.listA.txt,fusions.listB.txt \
                --genome_lib ${genome_lib} \
                --left_fq ${fastq[0]} --right_fq ${fastq[1]} \
                --out_dir . \
                --out_prefix finspector \
                --vis
"""
}
