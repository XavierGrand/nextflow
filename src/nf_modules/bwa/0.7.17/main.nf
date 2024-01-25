version = "0.7.17"
container_url = "lbmc/bwa:${version}"


workflow mapping {
  take:
    fasta
    fastq
  main:
    index_fasta(fasta)
    mapping_fastq(index_fasta.out.index.collect(), fastq)
  emit:
    bam = mapping_fastq.out.bam
    report = mapping_fastq.out.report
}


params.index_fasta = ""
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("${fasta.simpleName}.*"), emit: index
    tuple val(file_id), path("*_bwa_report.txt"), emit: report

  script:
"""
bwa index ${params.index_fasta} -p ${fasta.simpleName} ${fasta} \
&> ${fasta.simpleName}_bwa_report.txt
"""
}


params.mapping_fastq = ""
params.mapping_fastq_out = ""
params.file_suffix = ""
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.bam"), emit: bam
  tuple val(file_id), path("${file_prefix}_bwa_report.txt"), emit: report

  script:
  if (file_id instanceof List){
    library = file_id[0]
    file_prefix = file_id[0]
  } else if (file_id instanceof Map) {
      library = file_id[0]
      file_prefix = file_id[0]
      if (file_id.containsKey('library')) {
        library = file_id.library
        file_prefix = file_id.id
      }
  } else {
    library = file_id
    file_prefix = file_id
  }
bwa_mem_R = "@RG\\tID:${library}\\tSM:${library}\\tLB:lib_${library}\\tPL:illumina"
  if (reads.size() == 2)
"""
bwa mem -t ${task.cpus} \
${params.mapping_fastq} \
-R '${bwa_mem_R}' \
${index[0].baseName} ${reads[0]} ${reads[1]} 2> \
  ${file_prefix}_bwa_report.txt | \
  samtools view -@ ${task.cpus} -Sb - > ${file_prefix}${file_suffix}.bam
"""
  else
"""
bwa mem -t ${task.cpus} \
${params.mapping_fastq} \
-R '${bwa_mem_R}' \
${index[0].baseName} ${reads} 2> \
  ${file_prefix}_bwa_report.txt | \
  samtools view -@ ${task.cpus} -Sb - > ${file_prefix}${file_suffix}.bam
"""
}

