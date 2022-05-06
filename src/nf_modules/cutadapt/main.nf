version = "2.1"
container_url = "lbmc/cutadapt:${version}"

params.adapter_3_prim = "AGATCGGAAGAG"
params.adapter_5_prim = "CTCTTCCGATCT"
params.adaptor_removal = "-a ${params.adapter_3_prim} -g ${params.adapter_5_prim} -A ${params.adapter_3_prim} -G ${params.adapter_5_prim}"
params.adaptor_removal_out = ""
process adaptor_removal {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.adaptor_removal_out != "") {
    publishDir "results/${params.adaptor_removal_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*_cut_*"), emit: fastq
  path "*_report.txt", emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  if (reads.size() == 2)
  """
  cutadapt ${params.adaptor_removal} \
  -o ${file_prefix}_cut_R1.fastq.gz -p ${file_prefix}_cut_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${file_prefix}_report.txt
  """
  else
  """
  cutadapt ${params.adaptor_removal} \
  -o ${file_prefix}_cut.fastq.gz \
  ${reads} > ${file_prefix}_report.txt
  """
}

params.trim_quality = "20"
params.trimming = "-q ${params.trim_quality},${params.trim_quality}"
params.trimming_out = ""
process trimming {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.trimming_out != "") {
    publishDir "results/${params.trimming_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("*_trim_*"), emit:fastq
  path "*_report.txt", emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  if (reads.size() == 2)
  """
  cutadapt ${params.trimming} \
  -o ${file_prefix}_trim_R1.fastq.gz -p ${file_prefix}_trim_R2.fastq.gz \
  ${reads[0]} ${reads[1]} > ${file_prefix}_report.txt
  """
  else
  """
  cutadapt ${params.trimming} \
  -o ${file_prefix}_trim.fastq.gz \
  ${reads} > ${file_prefix}_report.txt
  """
}


process 5pRACE {
  container = "${container_url}"
  label "big_mem_mono_cpus"tag "$file_id"
  
  if (params.cutadapt_out != "") {
    publishDir "results/${params.cutadapt_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(fastq)

  output:
  tuple val(file_id), path("*_cut_*"), emit: fastq

  """
  cutadapt -e 0.2 -g CGACTGGAGCACGAGGACACTGACATGGACTGAAGGAGTAGAAA -g TTAGGCAGAGGTGAAAAAGTTG 
                  -a TTTCTACTCCTTCAGTCCATGTCAGTGTCCTCGTGCTCCAGTCG -a CAACTTTTTCACCTCTGCCTAA 
                  -o ${}
                  ${fastq}
  """
}