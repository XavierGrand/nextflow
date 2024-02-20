version = "2.25.0"
container_url = "lbmc/bedtools:${version}"

params.fasta_from_bed = "-name"
params.fasta_from_bed_out = ""
process fasta_from_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.fasta_from_bed_out != "") {
    publishDir "results/${params.fasta_from_bed_out}", mode: 'copy'
  }

  input:
  tuple val(fasta_id), path(fasta)
  tuple val(file_id), path(bed)

  output:
  tuple val(file_id), path("*_extracted.fasta"), emit: fasta

  script:
"""
bedtools getfasta ${params.fasta_from_bed} \
-fi ${fasta} -bed ${bed} -fo ${bed.baseName}_extracted.fasta
"""
}

params.merge_bed = ""
params.merge_bed_out = ""
process merge_bed {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.merge_bed_out != "") {
    publishDir "results/${params.merge_bed_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bed)

  output:
  tuple val(file_id), path("*_merged.fasta"), emit: bed

  script:
"""
bedtools merge ${params.merge_bed} -i ${bed} > ${bed[0].simpleName}_merged.bed
"""
}

params.bam_to_fastq_singleend = ""
params.bam_to_fastq_singleend_out = ""
process bam_to_fastq_singleend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"
  if (params.bam_to_fastq_singleend_out != "") {
    publishDir "results/${params.bam_to_fastq_singleend_out}", mode: 'copy'
  }

  input:
  tuple val(bam_id), path(bam)

  output:
  tuple val(bam_id), path("*.fastq"), emit: fastq

  script:
"""
bedtools bamtofastq \
  ${params.bam_to_fastq_singleend} \
  -i ${bam} -fq ${bam.baseName}.fastq
"""
}

params.bam_to_fastq_pairedend = ""
params.bam_to_fastq_pairedend_out = ""
process bam_to_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"
  if (params.bam_to_fastq_pairedend_out != "") {
    publishDir "results/${params.bam_to_fastq_pairedend_out}", mode: 'copy'
  }

  input:
  tuple val(bam_id), path(bam)

  output:
  tuple val(bam_id), path("*.fastq"), emit: fastq

  script:
"""
bedtools bamtofastq \
  ${params.bam_to_fastq_pairedend} \
  -i ${bam} -fq ${bam.baseName}_R1.fastq -fq2 ${bam.baseName}_R2.fastq
"""
}

params.bam_to_bedgraph = ""
params.bam_to_bedgraph_out = ""
process bam_to_bedgraph {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${bam_id}"
  if (params.bam_to_bedgraph_out != "") {
    publishDir "results/${params.bam_to_bedgraph_out}", mode: 'copy'
  }

  input:
  tuple val(bam_id), path(bam)

  output:
  tuple val(bam_id), path("*.bed"), emit: bed

  script:
"""
bedtools genomecov \
  ${params.bam_to_bedgraph} \
  -ibam ${bam} \
  -bg > ${bam.simpleName}.bed
"""
}
