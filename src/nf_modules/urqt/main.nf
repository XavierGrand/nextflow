version = "d62c1f8"
container_url = "lbmc/urqt:${version}"

trim_quality = "20"

params.trimming = "--t 20"
process trimming {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}"

  input:
  tuple val(file_id), path(reads)

  output:
  tuple val(pair_id), path("*_trim_R{1,2}.fastq.gz"), emit: fastq
  path "*_report.txt", emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  if (reads.size() == 2)
"""
UrQt ${params.trimming} --m ${task.cpus} --gz \
  --in ${reads[0]} --inpair ${reads[1]} \
  --out ${file_prefix}_trim_R1.fastq.gz --outpair ${file_prefix}_trim_R2.fastq.gz \
  > ${pair_id}_trimming_report.txt
"""
  else
"""
UrQt ${params.trimming} --m ${task.cpus} --gz \
  --in ${reads[0]} \
  --out ${file_prefix}_trim.fastq.gz \
  > ${file_prefix}_trimming_report.txt
"""
}