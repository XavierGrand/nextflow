version = "0.20.1"
container_url = "lbmc/fastp:${version}"

params.fastp_protocol = ""

params.fastp = ""
params.fastp_out = ""
workflow fastp {
  take:
    fastq

  main:
  switch(params.fastp_protocol) {
    case "accel_1splus":
      fastp_accel_1splus(fastq)
      fastp_accel_1splus.out.fastq.set{res_fastq}
      fastp_accel_1splus.out.report.set{res_report}
    break;
    default:
      fastp_default(fastq)
      fastp_default.out.fastq.set{res_fastq}
      fastp_default.out.report.set{res_report}
    break;
  }
  emit:
    fastq = res_fastq
    report = res_report
}

process fastp_default {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.fastp_out != "") {
    publishDir "results/${params.fastp_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
    tuple val(file_id), path("*_trim.fastq.gz"), emit: fastq
    tuple val(file_id), path("${file_prefix}.html"), emit: html
    tuple val(file_id), path("${file_prefix}_fastp.json"), emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  if (reads.size() == 2)
  """
  fastp --thread ${task.cpus} \
    --qualified_quality_phred 30 \
    --disable_length_filtering \
    --detect_adapter_for_pe \
    ${params.fastp} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]} \
    --out1 ${file_prefix}_R1_trim.fastq.gz \
    --out2 ${file_prefix}_R2_trim.fastq.gz \
    --html ${file_prefix}.html \
    --json ${file_prefix}_fastp.json \
    --report_title ${file_prefix}
  """
  else
  """
  fastp --thread ${task.cpus} \
    --qualified_quality_phred 20 \
    --disable_length_filtering \
    ${params.fastp} \
    --in1 ${reads[0]} \
    --out1 ${file_prefix}_trim.fastq.gz \
    --html ${file_prefix}.html \
    --json ${file_prefix}_fastp.json \
    --report_title ${file_prefix}
  """
}

process fastp_accel_1splus {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_prefix"
  if (params.fastp_out != "") {
    publishDir "results/${params.fastp_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
    tuple val(file_id), path("*_trim.fastq.gz"), emit: fastq
    tuple val(file_id), path("${file_prefix}.html"), emit: html
    tuple val(file_id), path("${file_prefix}_fastp.json"), emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }

  if (reads.size() == 2)
  """
  fastp --thread ${task.cpus} \
    --disable_quality_filtering \
    --disable_length_filtering \
    --disable_trim_poly_g \
    --detect_adapter_for_pe \
    --stdout \
    --in1 ${reads[0]} \
    --in2 ${reads[1]} 2> /dev/null | \
    fastp --thread ${task.cpus} \
      --stdin \
      --interleaved_in \
      --trim_front1=10 \
      --trim_front2=10 \
      --disable_adapter_trimming \
      --qualified_quality_phred 20 \
      --disable_length_filtering \
      --detect_adapter_for_pe \
      ${params.fastp} \
      --out1 ${file_prefix}_R1_trim.fastq.gz \
      --out2 ${file_prefix}_R2_trim.fastq.gz \
      --html ${file_prefix}.html \
      --json ${file_prefix}_fastp.json \
      --report_title ${file_prefix}
  """
  else
  """
  fastp --thread ${task.cpus} \
    --disable_quality_filtering \
    --disable_length_filtering \
    --disable_trim_poly_g \
    --detect_adapter_for_pe \
    --stdout \
    --in1 ${reads[0]} 2> /dev/null  | \
    fastp --thread ${task.cpus} \
      --disable_adapter_trimming \
      --stdin \
      --trim_front1=10 \
      --qualified_quality_phred 20 \
      --disable_length_filtering \
      --detect_adapter_for_pe \
      ${params.fastp} \
      --out1 ${file_prefix}_trim.fastq.gz \
      --html ${file_prefix}.html \
      --json ${file_prefix}_fastp.json \
      --report_title ${file_prefix}
  """
}

