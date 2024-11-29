version = "9.0.0"
container_url = "xgrand/cellranger:${version}"


params.cellranger_count_out = ""
process cellranger_count {
    tag ""
    label "huge_mem_multi_cpus"

    if (params.cellranger_count_out != "") {
        publishDir "results/${params.cellranger_count_out}", mode: 'copy'
    }

    input:
        tuple val(sample), path(fastq)
        path(refdata)

    output:
        path "merged.fastq.gz", emit: merged_fastq

    script:
    memory = "${task.memory}" - ~/\s*GB/
"""
cellranger count --id=${sample} \
	    --transcriptome=${refdata} \
	    --fastqs=${fastq} \
	    --sample=${sample} \
	    --localcores=${task.cpus} \
	    --localmem=${memory} \
	    --create-bam true
"""
}