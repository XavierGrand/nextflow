//Need to properly define a container
version = "1.2.3"
container_url = "lbmc/seekSV:${version}"

//Get soft-clipped reads from original bam file
params.get_soft_clipped = ""
params.soft_clipped_out = ""
params.soft_clipped_cache = ""
process get_soft_clipped {
  if (params.soft_clipped_cache != "") {
    cache false
  }
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.soft_clipped_out != "") {
    publishDir "results/${params.soft_clipped_out}", mode: 'copy'
  }

  input:  
    tuple val(file_id), path(bam)
    //tuple val(file_id), path(bam), path(bai)

  output:
    tuple val(file_id), path("*.clip.fq.gz"), emit: clip_fq
    tuple val(file_id), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_fq

  script:
  bam_only=bam[0][0]
  """
  seeksv getclip ${params.get_soft_clipped} -o ${file_id}_seeksv ${bam_only}
  """
}

/*
params.get_sv = ""
process get_sv {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:  
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.fq.gz*"), emit: clip_fq

  script:
"""
seeksv getsv ${params.get_sv} -o ${bam.baseName}_seeksv ${bam}
seeksv getsv /path/to/prefix.clip.bam \
             /path/to/input.bam \
             /path/to/prefix.clip.gz \
             /path/to/outputs/output.sv.txt \
             /path/to/outputs/output.unmapped.clip.fq.gz
"""
}
*/
