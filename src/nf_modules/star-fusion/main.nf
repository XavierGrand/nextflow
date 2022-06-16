version = "1.11.0"
// container_url = "lbmc/star-fusion:${version}"
container_url = "xgrand/star-fusion:${version}"

process star_fusion {

    container = "${container_url}"
    label "big_mem_multi_cpus"
    tag "$file_id"
    if (params.star-fusion_out != "") {
        publishDir "results/${params.star-fusion_out}", mode: 'copy'
    }

    input:
        path(genome_lib)
        tuple val(file_id), path(fastq)

    output:
        tuple val(), path ("star-fusion.fusion_predictions.tsv"), emit: fusion_predictions

    script:

if (fastq.size() == 2)
"""
STAR-Fusion --genome_lib_dir ${genome_lib} \
             --left_fq ${fastq[0]} \
             --right_fq ${fastq[1]} \
             --output_dir .
"""
else
"""
STAR-Fusion --genome_lib_dir ${genome_lib} \
             --left_fq ${fastq} \
             --output_dir .
"""
}
