version = "1.0"
container_url = "xgrand/viral-track:${version}"


params.viral-track_out = ""
process viral-track {
    tag "$sample"
    label "big_mem_multi_cpus"

    if (params.viral-track_out != "") {
        publishDir "results/${params.viral-track_out}", mode: 'copy'
    }

    input:
        val(sample) 
        path(fastq)
        path(refdata)

    output:
        path "", emit: 

    script:
    memory = "${task.memory}" - ~/\s*GB/
"""

"""
}