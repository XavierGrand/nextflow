version = "0.6.7"
container_url = "lbmc/sambamba:${version}"

params.index_bam = ""
process index_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
sambamba index ${params.index_bam} -t ${task.cpus} ${bam}
"""
}

params.sort_bam = ""
process sort_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
sambamba sort -t ${task.cpus} ${params.sort_bam} -o ${bam.baseName}_sorted.bam ${bam}
"""
}

params.mark_dup = ""
process mark_dup {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
sambamba markdup -t ${task.cpus} ${params.mark_dup} -o ${bam.baseName}_markdup.bam ${bam}
"""
}


params.split_bam = ""
process split_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_forward.bam*"), emit: bam_forward
    tuple val(file_id), path("*_reverse.bam*"), emit: bam_reverse
  script:
"""
sambamba view -t ${task.cpus} ${params.split_bam} -h -F "strand == '+'" ${bam} > \
  ${bam.baseName}_forward.bam
sambamba view -t ${task.cpus} ${params.split_bam} -h -F "strand == '-'" ${bam} > \
  ${bam.baseName}_reverse.bam
"""
}
