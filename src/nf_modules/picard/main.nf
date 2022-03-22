version = "2.18.11"
container_url = "lbmc/picard:${version}"

params.mark_duplicate = "VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true"
params.mark_duplicate_out = ""
process mark_duplicate {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.mark_duplicate_out != "") {
    publishDir "results/${params.mark_duplicate_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)
  output:
    tuple val(file_id) , path("*.bam"), emit: bam
    path "*_report.dupinfo.txt", emit: report


  script:
"""
PicardCommandLine MarkDuplicates \
  ${params.mark_duplicate} \
  INPUT=${bam} \
  OUTPUT=${bam.baseName}_dedup.bam \
  METRICS_FILE=${bam.baseName}_picard_dedup_report.dupinfo.txt &> \
  picard_${bam.baseName}.log
"""
}

params.normalize_fasta = ""
params.normalize_fasta_out = ""
process normalize_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.normalize_fasta_out != "") {
    publishDir "results/${params.normalize_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
  output:
    tuple val(file_id), path("results/*.fasta.gz"), emit: fasta 

  script:
"""
mkdir -p results
PicardCommandLine NormalizeFasta \
      I=${fasta} \
      O=results/${fasta.simpleName}.fasta
gzip results/${fasta.simpleName}.fasta
"""
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
    tuple val(file_id), path("*.dict"), emit: index

  script:
"""
PicardCommandLine CreateSequenceDictionary \
  ${params.index_fasta} \
  REFERENCE=${fasta} \
  OUTPUT=${fasta.baseName}.dict
"""
}

params.index_bam = ""
params.index_bam_out = ""
process index_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_bam_out != "") {
    publishDir "results/${params.index_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)
  output:
    tuple val(file_id), path("*"), emit: index

  script:
"""
PicardCommandLine BuildBamIndex \
  ${params.index_bam} \
  INPUT=${bam}
"""
}
