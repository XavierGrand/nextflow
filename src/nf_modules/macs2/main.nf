version = "2.1.2"
container_url = "lbmc/macs2:${version}"

params.macs_gsize=3e9
params.macs_mfold="5 50"
params.peak_calling = "--mfold ${params.macs_mfold} --gsize ${params.macs_gsize}"
params.peak_calling_out = ""
process peak_calling {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.peak_calling_out != "") {
    publishDir "results/${params.peak_calling_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam_ip), path(bam_control)

  output:
    tuple val(file_id), path("*.narrowPeak"), emit: peak
    tuple val(file_id), path("*.bed"), emit: summits
    tuple val(file_id), path("*_peaks.xls"), path("*_report.txt"), emit: report

  script:
/* remove --nomodel option for real dataset */
"""
macs2 callpeak \
  ${params.peak_calling} \
  --treatment ${bam_ip} \
  --call-summits \
  --control ${bam_control} \
  --keep-dup all \
  --qvalue 0.99 \
  --name ${bam_ip.simpleName} 2> \
  ${bam_ip.simpleName}_macs2_report.txt

if grep -q "ERROR" ${bam_ip.simpleName}_macs2_report.txt; then
  echo "MACS3 error"
  exit 1
fi
"""
}

params.peak_calling_bg = "--mfold ${params.macs_mfold} --gsize ${params.macs_gsize}"
params.peak_calling_bg_out = ""
process peak_calling_bg {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.peak_calling_bg_out != "") {
    publishDir "results/${params.peak_calling_bg_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bg_ip), path(bg_control)

  output:
    tuple val(file_id), path("*.narrowPeak"), emit: peak
    tuple val(file_id), path("*.bed"), emit: summits
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
/* remove --nomodel option for real dataset */
"""
awk '{print \$1"\t"\$2"\t"\$3"\t.\t+\t"\$4}' ${bg_ip} > \
  ${bg_ip.simpleName}.bed
awk '{print \$1"\t"\$2"\t"\$3"\t.\t+\t"\$4}' ${bg_control} > \
  ${bg_control.simpleName}.bed
macs2 callpeak \
  ${params.peak_calling_bg} \
  --treatment ${bg_ip.simpleName}.bed \
  --qvalue 0.99 \
  --call-summits \
  --control ${bg_control.simpleName}.bed \
  --keep-dup all \
  --name ${bg_ip.simpleName} 2> \
  ${bg_ip.simpleName}_macs2_report.txt

if grep -q "ERROR" ${bg_ip.simpleName}_macs2_report.txt; then
  echo "MACS3 error"
  exit 1
fi
"""
}

