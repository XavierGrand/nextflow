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
    tuple val(file_id), path("*.clip.gz"), emit: clip_gz
    tuple val(file_id), path("*.unmapped_{1,2}.fq.gz"), emit: unmapped_fq

  script:
  bam_only=bam[0][0]
  """
  seeksv getclip ${params.get_soft_clipped} -o ${file_id}_seeksv ${bam_only}
  """
}

params.get_sv = ""
params.get_sv_out = ""
process get_sv {
  if (params.get_sv_cache != "") {
    cache false
  }
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  if (params.get_sv_out != "") {
    publishDir "results/${params.get_sv_out}", mode: 'copy'
  }

  input:  
    tuple val(file_id), path(clip_bam)
    tuple val(file_id), path(original_bam)
    tuple val(file_id), path(clip_gz)


  output:
    tuple val(file_id), path("*.sv.txt"), emit: sv_report
    tuple val(file_id), path("*.clip.fq.gz"), emit: unmapped_fq

  script:
  if (original_bam instanceof List ){
    original_bam_file = original_bam[0]
  } else {
    original_bam_file = original_bam
  }
  """
  seeksv getsv ${params.get_sv} ${clip_bam} ${original_bam_file} ${clip_gz} ${file_id}/${file_id}_seekSV.sv.txt ${file_id}/${file_id}_seekSV.unmapped.clip.fq.gz  
  """
}
