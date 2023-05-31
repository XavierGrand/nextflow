version = "3.15.4"
container_url = "xgrand/alpine:${version}"


params.fastq_out = ""
process concatenate {
    tag "fastq_folder"
    label "big_mem_mono_cpus"

    if (params.fastq_out != "") {
        publishDir "results/${params.fastq_out}", mode: 'copy'
    }

    input:
        path fastq

    output:
        path "merged.fastq.gz", emit: merged_fastq

    script:
"""
zcat ${fastq}/*.fastq.gz > merged.fastq
gzip --quiet merged.fastq
"""
}

params.concat_fusion_out = ""
process concat_fusion {
    tag "concat_fusion"
    label "big_mem_mono_cpus"

    if (params.fastq_out != "") {
        publishDir "results/${concat_fusion_out}", mode: 'copy'
    }

    input:
        path fusions
        path discarded_fusions

    output:
        path "*_fusions.tsv", emit: concatenated_fusions

    script:
"""
cat ${fusions} > res_fusions.tsv
tail +2 ${discarded_fusions} >> res_fusions.tsv
"""
}