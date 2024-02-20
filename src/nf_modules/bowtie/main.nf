version = "1.2.2"
container_url = "lbmc/bowtie:${version}"

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
    tuple val(file_id), path("*.index*"), emit: index
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
bowtie-build --threads ${task.cpus} \
  ${params.index_fasta} \
  -f ${fasta} ${fasta.baseName}.index &> \
  ${fasta.baseName}_bowtie_index_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie_index_report.txt; then
  exit 1
fi
"""
}

params.mapping_fastq = "--very-sensitive"
params.mapping_fastq_out = ""
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
  if (reads.size() == 2)
  """
  # -v specify the max number of missmatch, -k the number of match reported per
  # reads
  bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
    ${params.mapping_fastq} \
    -1 ${reads[0]} -2 ${reads[1]} 2> \
    ${file_id}_bowtie_report_tmp.txt | \
    samtools view -Sb - > ${file_id}.bam

  if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
    exit 1
  fi
  tail -n 19 ${file_id}_bowtie_report_tmp.txt > \
    ${file_id}_bowtie_mapping_report.txt
  """
  else
  """
  bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
    ${params.mapping_fastq}
    -q ${reads} 2> \
    ${file_id}_bowtie_report_tmp.txt | \
    samtools view -Sb - > ${file_id}.bam

  if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
    exit 1
  fi
  tail -n 19 ${file_id}_bowtie_report_tmp.txt > \
    ${file_id}_bowtie_mapping_report.txt
  """
}

params.mapping_fastq_pairedend = ""
process mapping_fastq_pairedend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$pair_id"

  input:
  path index
  tuple val(pair_id), path(reads)

  output:
  tuple val(pair_id), path("*.bam"), emit: bam
  path "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
# -v specify the max number of missmatch, -k the number of match reported per
# reads
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
  ${params.mapping_fastq_pairedend} \
  -1 ${reads[0]} -2 ${reads[1]} 2> \
  ${pair_id}_bowtie_report_tmp.txt | \
  samtools view -Sb - > ${pair_id}.bam

if grep -q "Error" ${pair_id}_bowtie_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${pair_id}_bowtie_report_tmp.txt > \
  ${pair_id}_bowtie_mapping_report.txt
"""
}

params.mapping_fastq_singleend = ""
process mapping_fastq_singleend {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
  path index
  tuple val(file_id), path(reads)

  output:
  set file_id, "*.bam", emit: bam
  file "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
        index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
    }
  }
"""
bowtie --best -v 3 -k 1 --sam -p ${task.cpus} ${index_id} \
  ${params.mapping_fastq_singleend} \
  -q ${reads} 2> \
  ${file_id}_bowtie_report_tmp.txt | \
  samtools view -Sb - > ${file_id}.bam

if grep -q "Error" ${file_id}_bowtie_report_tmp.txt; then
  exit 1
fi
tail -n 19 ${file_id}_bowtie_report_tmp.txt > \
  ${file_id}_bowtie_mapping_report.txt
"""
}
