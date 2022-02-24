version = "0.26.0"
container_url = "lbmc/kb:${version}"

params.index_fasta = ""
params.index_fasta_out = ""

workflow index_fasta {
  take:
    fasta
    gtf

  main:
    tr2g(gtf)
    index_default(fasta, gtf, tr2g.out.t2g)

  emit:
    index = index_default.out.index
    t2g = index_default.out.t2g
    report = index_default.out.report
}

process tr2g {
  // create transcript to gene table from gtf if no transcript to gene file is provided
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gtf)

  output:
    tuple val(file_id), path("t2g.txt"), emit: t2g

  script:
  """
  t2g.py --gtf ${gtf}
  sort -k1 -u t2g_dup.txt > t2g.txt
  """
}

process g2tr {
  // create gene to transcript table from gtf if no transcript to gene file is provided
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gtf)

  output:
    tuple val(file_id), path("g2t.txt"), emit: g2t

  script:
  """
  t2g.py --gtf ${gtf}
  sort -k1 -u t2g_dup.txt > t2g.txt
  awk 'BEGIN{OFS="\\t"}{print \$2, \$1}' t2g.txt > g2t.txt
  """
}

process index_default {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
    tuple val(gtf_id), path(gtf)
    tuple val(t2g_id), path(transcript_to_gene)

  output:
    tuple val(file_id), path("*.idx"), emit: index
    tuple val(t2g_id), path("${transcript_to_gene}"), emit: t2g
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
kb ref \
  -i ${fasta.simpleName}.idx \
  -g ${transcript_to_gene} \
  ${params.index_fasta} \
  -f1 cdna.fa ${fasta} ${gtf} > ${fasta.simpleName}_kb_index_report.txt
"""
}


include { split } from "./../flexi_splitter/main.nf"

params.kb_protocol = "10x_v3"
params.count = ""
params.count_out = ""
workflow count {
  take:
    index
    fastq
    transcript_to_gene
    whitelist
    config

  main:
  whitelist
    .ifEmpty(["NO WHITELIST", 0])
    .set{ whitelist_optional }
  switch(params.kb_protocol) {
    case "marsseq":
      split(fastq, config.collect())
      kb_marseq(index.collect(), split.out.fastq, transcript_to_gene.collect(), whitelist_optional.collect())
      kb_marseq.out.counts.set{res_counts}
      kb_marseq.out.report.set{res_report}
    break;
    default:
      kb_default(index.collect(), fastq, transcript_to_gene.collect(), whitelist_optional.collect())
      kb_default.out.counts.set{res_counts}
      kb_default.out.report.set{res_report}
    break;
  }

  emit:
    counts = res_counts
    report = res_report
}

process kb_default {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.count_out != "") {
    publishDir "results/${params.count_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)
  tuple val(t2g_id), path(transcript_to_gene)
  tuple val(whitelist_id), path(whitelist)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  def kb_memory = "${task.memory}" - ~/GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO WHITELIST"){
    whitelist_param = "-w ${whitelist}"
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${kb_memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    -o ${file_prefix} \
    ${whitelist_param} \
    -x 10XV3 \
    --h5ad \
    ${params.count} \
    ${reads[0]} ${reads[1]} > ${file_prefix}_kb_mapping_report.txt
  
  fix_t2g.py --t2g ${transcript_to_gene}
  cp fix_t2g.txt ${file_prefix}/
  cp ${transcript_to_gene} ${file_prefix}/
  """
}

process kb_marseq {
  // With the MARS-Seq protocol, we have:
  // on the read 1: 4 nt of bc plate
  // on the read 2: 6 nt of bc cell, and 8 nt of UMI
  // this process expect that the bc plate is removed from the read 1
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.count_out != "") {
    publishDir "results/${params.count_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)
  tuple val(t2g_id), path(transcript_to_gene)
  tuple val(whitelist_id), path(whitelist)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  def kb_memory = "${task.memory}" - ~/GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO WHITELIST"){
    whitelist_param = "-w ${whitelist}"
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${kb_memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    -o ${file_prefix} \
    ${whitelist_param} \
    ${params.count} \
    --h5ad \
    -x 1,0,6:1,6,14:0,0,0 \
    ${reads[0]} ${reads[1]} > ${file_prefix}_kb_mapping_report.txt
  fix_t2g.py --t2g ${transcript_to_gene}
  cp fix_t2g.txt ${file_prefix}/
  cp ${transcript_to_gene} ${file_prefix}/
  """
  else
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${kb_memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    -o ${file_prefix} \
    ${whitelist_param} \
    ${params.count} \
    -x 1,0,6:1,6,14:0,0,0 \
    --h5ad \
    ${reads} > ${file_prefix}_kb_mapping_report.txt
  fix_t2g.py --t2g ${transcript_to_gene}
  cp fix_t2g.txt ${file_prefix}/
  cp ${transcript_to_gene} ${file_prefix}/
  """
}

// ************************** velocity workflow **************************

workflow index_fasta_velocity {
  take:
    fasta
    gtf

  main:
    tr2g(gtf)
    index_fasta_velocity_default(fasta, gtf, tr2g.out.t2g)

  emit:
    index = index_fasta_velocity_default.out.index
    t2g = index_fasta_velocity_default.out.t2g
    report = index_fasta_velocity_default.out.report
}

process index_fasta_velocity_default {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
    tuple val(gtf_id), path(gtf)
    tuple val(t2g_id), path(transcript_to_gene)

  output:
    tuple val(file_id), path("*.idx"), emit: index
    tuple val(t2g_id), path("${transcript_to_gene}"), path("cdna_t2c.txt"), path("intron_t2c.txt"), emit: t2g
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
kb ref \
  -i ${fasta.simpleName}.idx \
  -g ${transcript_to_gene} \
  ${params.index_fasta} \
  -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno \
  ${fasta} ${gtf} > ${fasta.simpleName}_kb_index_report.txt
"""
}

params.count_velocity = ""
params.count_velocity_out = ""
workflow count_velocity {
  take:
    index
    fastq
    transcript_to_gene
    whitelist
    config

  main:
  whitelist
    .ifEmpty(["NO WHITELIST", 0])
    .set{ whitelist_optional }
  switch(params.kb_protocol) {
    case "marsseq":
      split(fastq, config.collect())
      velocity_marseq(index.collect(), split.out.fastq, transcript_to_gene.collect(), whitelist_optional.collect())
      velocity_marseq.out.counts.set{res_counts}
      velocity_marseq.out.report.set{res_report}
    break;
    default:
      velocity_default(index.collect(), fastq, transcript_to_gene.collect(), whitelist_optional.collect())
      velocity_default.out.counts.set{res_counts}
      velocity_default.out.report.set{res_report}
    break;
  }

  emit:
    counts = res_counts
    report = res_report
}

process velocity_default {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.count_velocity_out != "") {
    publishDir "results/${params.count_velocity_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)
  tuple val(t2g_id), path(transcript_to_gene), path(cdna_t2g), path(intron_t2g)
  tuple val(whitelist_id), path(whitelist)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  def kb_memory = "${task.memory}" - ~/GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO WHITELIST"){
    whitelist_param = "-w ${whitelist}"
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${kb_memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    -o ${file_prefix} \
    -c1 ${cdna_t2g} \
    -c2 ${intron_t2g} \
    --workflow lamanno \
    ${whitelist_param} \
    -x 10XV3 \
    --h5ad \
    ${params.count} \
    ${reads[0]} ${reads[1]} > ${file_prefix}_kb_mapping_report.txt
  fix_t2g.py --t2g ${transcript_to_gene}
  cp fix_t2g.txt ${file_prefix}/
  cp ${transcript_to_gene} ${file_prefix}/
  cp ${cdna_t2g} ${file_prefix}/
  cp ${intron_t2g} ${file_prefix}/
  """
}

process velocity_marseq {
  // With the MARS-Seq protocol, we have:
  // on the read 1: 4 nt of bc plate
  // on the read 2: 6 nt of bc cell, and 8 nt of UMI
  // this process expect that the bc plate is removed from the read 1
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.count_velocity_out != "") {
    publishDir "results/${params.count_velocity_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)
  tuple val(t2g_id), path(transcript_to_gene), path(cdna_t2g), path(intron_t2g)
  tuple val(whitelist_id), path(whitelist)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  def kb_memory = "${task.memory}" - ~/GB/
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  def whitelist_param = ""
  if (whitelist_id != "NO WHITELIST"){
    whitelist_param = "-w ${whitelist}"
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${kb_memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    -o ${file_prefix} \
    -c1 ${cdna_t2g} \
    -c2 ${intron_t2g} \
    --workflow lamanno \
     --h5ad \
    ${whitelist_param} \
    ${params.count} \
    -x 1,0,6:1,6,14:0,0,0 \
    ${reads[0]} ${reads[1]} > ${file_prefix}_kb_mapping_report.txt
  fix_t2g.py --t2g ${transcript_to_gene}
  cp fix_t2g.txt ${file_prefix}/
  cp ${transcript_to_gene} ${file_prefix}/
  cp ${cdna_t2g} ${file_prefix}/
  cp ${intron_t2g} ${file_prefix}/
  """
  else
  """
  mkdir ${file_prefix}
  kb count  -t ${task.cpus} \
    -m ${kb_memory} \
    -i ${index} \
    -g ${transcript_to_gene} \
    -o ${file_prefix} \
    -c1 ${cdna_t2g} \
    -c2 ${intron_t2g} \
    --workflow lamanno \
    ${whitelist_param} \
    ${params.count} \
    -x 1,0,6:1,6,14:0,0,0 \
    ${reads} > ${file_prefix}_kb_mapping_report.txt
  fix_t2g.py --t2g ${transcript_to_gene}
  cp fix_t2g.txt ${file_prefix}/
  cp ${transcript_to_gene} ${file_prefix}/
  cp ${cdna_t2g} ${file_prefix}/
  cp ${intron_t2g} ${file_prefix}/
  """
}

