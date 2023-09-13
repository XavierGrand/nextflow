version = "2.2.1"
container_url = "lbmc/hisat2:${version}"

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
    tuple val(file_id), path("*.ht2*"), emit: index
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
  filename = fasta.toString()
  extention = filename[filename.lastIndexOf('.')..filename.size() - 1]
  if ( extention == ".gz" ) {
    base_name = fasta.baseName
    simple_name = fasta.simpleName
  } else {
    base_name = filename
    simple_name = fasta.simpleName
  }
"""
if [[ ${filename} == *.gz ]]; then
  gunzip ${filename}
fi
hisat2-build -p ${task.cpus} \
  ${base_name} \
  ${simple_name} &> \
  ${simple_name}_hisat2_index_report.txt

if grep -q "Error" ${simple_name}_hisat2_index_report.txt; then
  exit 1
fi
"""
}

params.mapping_fastq = ""
params.mapping_fastq_out = ""
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
  path "*_report.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.ht2.*/) {
        index_id = ( index_file =~ /(.*)\.1\.ht2.*/)[0][1]
    }
  }
  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  if (reads.size() == 2)
  """
  hisat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -1 ${reads[0]} \
    -2 ${reads[1]} 2> \
    ${file_prefix}_ht2_mapping_report.txt \
    | samtools view -@ ${task.cpus} -bS - \
    | samtools sort -@ ${task.cpus} -o ${file_prefix}.bam

  if grep -q "Error" ${file_prefix}_ht2_mapping_report.txt; then
    exit 1
  fi
  """
  else
  """
  hisat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -U ${reads} 2> \
    ${file_prefix}_ht2_mapping_report.txt \
    | samtools view -@ ${task.cpus} -bS - \
    | samtools sort -@ ${task.cpus} -o ${file_prefix}.bam
  if grep -q "Error" ${file_prefix}_ht2_mapping_report.txt; then
    exit 1
  fi
  """
}


params.hisat2 = ""
params.notaligned_name = "notaligned"
process genome_mapping {
  tag "$file_id"
  label "big_mem_multi_cpus"
  container = "${container_url}"
  publishDir "results/${params.folder}", mode: 'copy'

  input:
    tuple val(file_id), file(fastq_filtred)
    tuple val(index_id), path(index)

  output:
    tuple val(file_id), path("*${params.notaligned_name}*.fastq.gz"), emit: unaligned
    tuple val(file_id), path("*_aligned*.bam"), emit: aligned
    path "*.txt", emit: report

  script:
    if (file_id instanceof List){
      file_prefix = file_id[0]
    } else {
      file_prefix = file_id
    }
    index_id = index[0]
    for (index_file in index) {
        if (index_file =~ /.*\.1\.ht2/ && !(index_file =~ /.*\.rev\.1\.ht2/)) {
            index_id = ( index_file =~ /(.*)\.1\.ht2/)[0][1]
        }
    }
    if (fastq_filtred.size() == 2)
      """
      hisat2 -x ${index_id} -p ${task.cpus} \
      -1 ${fastq_filtred[0]} -2  ${fastq_filtred[1]} \
      --un-conc-gz ${file_prefix}_notaligned.fastq.gz \
      ${params.hisat2} \
      2> ${file_prefix}_hisat2_report.txt | samtools view  -@ ${task.cpus} -bS -f 2 -F 268 -q 10 -o ${file_prefix}_aligned.bam

      mv ${file_prefix}_notaligned.fastq.1.gz ${file_prefix}_${params.notaligned_name}_1.fastq.gz
      mv ${file_prefix}_notaligned.fastq.2.gz ${file_prefix}_${params.notaligned_name}_2.fastq.gz

      if grep -q "Error" ${file_prefix}_hisat2_report.txt; then
      exit 1
      fi
      """
    else
      """
      hisat2 -x ${index_id} -p ${task.cpus} \
      -U ${fastq_filtred} --un-gz ${file_prefix}_${params.notaligned_name}.fastq.gz \
      \
      2> ${file_prefix}_hisat2_report.txt | samtools view  -@ ${task.cpus} -bS -F 260 -q 10 -o ${file_prefix}_aligned.bam

      if grep -q "Error" ${file_prefix}_hisat2_report.txt; then
      exit 1
      fi
      """
}