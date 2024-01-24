version = "2.17"
container_url = "lbmc/minimap2:${version}"

params.index_fasta = ""
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("${fasta}"), path("*.mmi*"), emit: index

  script:
  memory = "${task.memory}" - ~/\s*GB/
"""
minimap2 ${params.index_fasta} -t ${task.cpus} -I ${memory}G -d ${fasta.baseName}.mmi ${fasta}
"""
}

params.mapping_fastq = "-ax sr"
params.mapping_fastq_out = ""
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(fasta_id), path(fasta), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.bam"), emit: bam

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  memory = "${task.memory}" - ~/\s*GB/
  memory = memory.toInteger() / (task.cpus + 1.0)
  if (reads.size() == 2)
  """
  minimap2 ${params.mapping_fastq} -t ${task.cpus} -K ${memory} ${fasta} ${reads[0]} ${reads[1]} |
    samtools view -Sb - > ${pair_id}.bam
  """
  else
  """
  minimap2 ${params.mapping_fastq} -t ${task.cpus} -K ${memory} ${fasta} ${reads} |
    samtools view -Sb - > ${pair_id}.bam
  """
}

params.mapping_hbv_genome = "-ax map-ont"
process mapping_hbv_genome {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${barcode}"
  if (params.mapping_hbv_genome_out != "") {
    publishDir "results/${params.mapping_hbv_genome_out}", mode: 'copy'
  }

  input:
  tuple val(barcode), path(fastq)
  tuple val(ref), path(genome)

  output:
  tuple val(barcode), path("${barcode}/${barcode}_res.bam"), emit: bam

  script:
  memory = "${task.memory}" - ~/\s*GB/
  memory = memory.toInteger() / (task.cpus + 1.0)
  """
  mkdir ${barcode}
  cd ${barcode}/
  minimap2 ${params.mapping_hbv_genome} -t ${task.cpus} -K ${memory} ../${genome} ../${fastq} |
    samtools view -Shb -F4 - > ${barcode}_res.bam
  """
}