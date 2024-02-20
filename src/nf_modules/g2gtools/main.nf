version = "0.2.8"
container_url = "lbmc/g2gtools:${version}"

params.vci_build = ""
params.vci_build_out = ""
process vci_build {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.vci_build_out != "") {
    publishDir "results/${params.vci_build_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta)
  output:
    tuple val(file_id), path("*.vci.gz"), path("*.vci.gz.tbi"), emit: vci
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }

  input_vcf = ""
  for (vcf_file in vcf) {
    input_vcf += " -i ${vcf_file}"
  }
"""
g2gtools vcf2vci \
  ${params.vci_build} \
  -p ${task.cpus} \
  -f ${fasta} \
  ${input_vcf} \
  -s ${file_prefix} \
  -o ${file_prefix}.vci 2> ${file_prefix}_g2gtools_vcf2vci_report.txt
"""
}

params.incorporate_snp = ""
params.incorporate_snp_out = ""
process incorporate_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.incorporate_snp_out != "") {
    publishDir "results/${params.incorporate_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vci), path(tbi)
    tuple val(ref_id), path(fasta)
  output:
    tuple val(file_id), path("${file_prefix}_snp.fa"), path("${vci}"), path("${tbi}"), emit: fasta
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
g2gtools patch \
  ${params.incorporate_snp} \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_prefix}_snp.fa 2> ${file_prefix}_g2gtools_path_report.txt
"""
}

params.incorporate_indel = ""
params.incorporate_indel_out = ""
process incorporate_indel {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.incorporate_indel_out != "") {
    publishDir "results/${params.incorporate_indel_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta), path(vci), path(tbi)
  output:
    tuple val(file_id), path("${file_prefix}_snp_indel.fa"), path("${vci}"), path("${tbi}"), emit: fasta
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
g2gtools transform \
  ${params.incorporate_indel} \
  -p ${task.cpus} \
  -i ${fasta} \
  -c ${vci} \
  -o ${file_prefix}_snp_indel.fa 2> ${file_prefix}_g2gtools_transform_report.txt
"""
}

params.convert_gtf = ""
params.convert_gtf_out = ""
process convert_gtf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.convert_gtf_out != "") {
    publishDir "results/${params.convert_gtf_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vci), path(tbi)
    tuple val(annot_id), path(gtf)
  output:
    tuple val(file_id), path("${file_prefix}.gtf"), emit: gtf
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
g2gtools convert \
  ${params.convert_gtf} \
  -i ${gtf} \
  -c ${vci} \
  -o ${file_prefix}.gtf 2> ${file_prefix}_g2gtools_convert_report.txt
"""
}

params.convert_bed = ""
params.convert_bed_out = ""
process convert_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.convert_bed_out != "") {
    publishDir "results/${params.convert_bed_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vci), path(tbi)
    tuple val(annot_id), path(bed)
  output:
    tuple val(file_id), path("${file_id}.bed"), emit: bed
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
g2gtools convert \
  ${params.convert_bed} \
  -i ${bed} \
  -c ${vci} \
  -o ${file_id}.bed 2> ${file_id}_g2gtools_convert_report.txt
"""
}

params.convert_bam = ""
params.convert_bam_out = ""
process convert_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id} ${file_id}"
  if (params.convert_bam_out != "") {
    publishDir "results/${params.convert_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vci), path(tbi)
    tuple val(bam_id), path(bam)
  output:
    tuple val(file_id), path("${file_id}_${bam_id.baseName}.bam"), emit: bam
    tuple val(file_id), path("*_report.txt"), emit: report
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
g2gtools convert \
  ${params.convert_bam} \
  -i ${bam} \
  -c ${vci} \
  -o ${file_id}_${bam.baseName}.bam 2> ${file_id}_g2gtools_convert_report.txt
"""
}