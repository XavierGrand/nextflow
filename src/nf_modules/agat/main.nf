version = "0.8.0"
container_url = "lbmc/agat:${version}"

params.gff_to_bed = ""
params.gff_to_bed_out = ""
process gff_to_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.gff_to_bed_out != "") {
    publishDir "results/${params.gff_to_bed_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gff)
  output:
    tuple val(file_id), path("*.bed"), emit: bed

  script:
"""
zcat ${gff} > ${gff.baseName}.gff
agat_convert_sp_gff2bed.pl ${params.gff_to_bed} --gff ${gff.baseName}.gff -o ${gff.simpleName}.bed
"""
}

params.gff_to_gtf = ""
params.gff_to_gtf_out = ""
process gff_to_gtf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.gff_to_gtf_out != "") {
    publishDir "results/${params.gff_to_gtf_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gff)
  output:
    tuple val(file_id), path("*.gtf"), emit: gtf

  script:
"""
zcat ${gff} > ${gff.baseName}.gff
agat_convert_sp_gff2gtf.pl ${params.gff_to_gtf} --gff ${gff.baseName}.gff -o ${gff.simpleName}.gtf
"""
}