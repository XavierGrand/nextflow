version = "3.5.1"
container_url = "lbmc/deeptools:${version}"

params.index_bam = ""
params.index_bam_out = ""
process index_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_bam_out != "") {
    publishDir "results/${params.index_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("${bam}"), path("*.bam*"), emit: bam_idx

  script:
"""
sambamba index -t ${task.cpus} ${bam}
"""
}

params.bam_to_bigwig = ""
params.bam_to_bigwig_out = ""
process bam_to_bigwig {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.bam_to_bigwig_out != "") {
    publishDir "results/${params.bam_to_bigwig_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam), path(idx)

  output:
    tuple val(file_id), path("*.bw"), emit: bw

  script:
"""
bamCoverage -p ${task.cpus} --ignoreDuplicates -b ${bam} \
  -o ${bam.simpleName}.bw
"""
}

params.compute_matrix = ""
params.compute_matrix_out = ""
process compute_matrix {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${bed_file_id}"
  if (params.compute_matrix_out != "") {
    publishDir "results/${params.compute_matrix_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bw)
    tuple val(bed_file_id), path(bed)

  output:
    tuple val(bed_file_id), path("*.mat.gz"), emit: matrix

  script:
"""
computeMatrix scale-regions -S ${bw} \
  -p ${task.cpus} \
  -R ${bed} \
  --beforeRegionStartLength 100 \
  --afterRegionStartLength 100 \
  -o ${bed.simpleName}.mat.gz
"""
}

params.plot_profile = ""
params.plot_profile_out = ""
process plot_profile {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.compute_matrix_out != "") {
    publishDir "results/${params.compute_matrix_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(matrix)

  output:
    tuple val(file_id), path("*.pdf"), emit: pdf

  script:
/*
see more option at
https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
*/
"""
plotProfile -m ${matrix} \
  --plotFileFormat=pdf \
  -out ${matrix.simpleName}.pdf \
  --plotType=fill \
  --perGroup \
  --plotTitle "${params.title}"
"""
}
