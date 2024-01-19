//Need to properly define a container
version = "1.2.1"
container_url = "lbmc/seekSV:${version}"

//Get soft-clipped reads from original bam file
params.get_soft_clipped = ""
params.soft_clipped_out = ""
process get_soft_clipped {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.soft_clipped_out != "") {
    publishDir "results/${params.soft_clipped_out}", mode: 'copy'
  }

  input:  
    //tuple val(file_id), path(bam)
    tuple val(file_id),path(bam),path(bai)

  output:
    tuple val(file_id), path("*.fq.gz*"), emit: clip_fq

  script:
//seeksv getclip ${params.get_soft_clipped} -o ${bam.baseName}_seeksv ${bam}
"""
seeksv getclip ${params.get_soft_clipped} -o ${bam.baseName}_seeksv ${bam}
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
