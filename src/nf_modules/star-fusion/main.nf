version = "1.11.0"
// container_url = "lbmc/star-fusion:${version}"
container_url = "xgrand/star-fusion:${version}"

params.star_fusion_out = "predicted_fusion"
process star_fusion {

    container = "${container_url}"
    label "mid_mem_multi_cpus"
    tag "$file_id"
    if (params.star_fusion_out != "") {
        publishDir "results/${params.star_fusion_out}", mode: 'copy'
    }

    input:
        path(genome_lib)
        tuple val(file_id), path(fastq)

    output:
        tuple val(file_id), path ("star-fusion.fusion_predictions.tsv"), emit: fusion_predictions

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
