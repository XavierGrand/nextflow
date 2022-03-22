version = "0.12.2"
container_url = "lbmc/gffread:${version}"

params.gffread = ""
params.gffread_out = ""
process gffread {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_prefix"
  if (params.gffread_out != "") {
    publishDir "results/${params.gffread_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(gtf)
  tuple val(fasta_id), path(fasta)

  output:
    tuple val(fasta_id), path("${file_prefix}.fasta"), emit: fasta

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  """
  gffread ${gtf} -g ${fasta} -M -x dup_${file_prefix}.fasta
  awk 'BEGIN {i = 1;} { if (\$1 ~ /^>/) { tmp = h[i]; h[i] = \$1; } else if (!a[\$1]) { s[i] = \$1; a[\$1] = "1"; i++; } else { h[i] = tmp; } } END { for (j = 1; j < i; j++) { print h[j]; print s[j]; } }' < dup_${file_prefix}.fasta | grep -v -e "^\$" > ${file_prefix}.fasta
  """
}

params.spliced_cds = ""
params.spliced_cds_out = ""
process spliced_cds {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_prefix"
  if (params.spliced_cds_out != "") {
    publishDir "results/${params.spliced_cds_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(gtf)
  tuple val(fasta_id), path(fasta)

  output:
    tuple val(fasta_id), path("${file_prefix}.fasta"), emit: fasta

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  """
gzip -dc ${fasta} > ${fasta.simpleName}_un.fasta
gzip -dc ${gtf} > ${gff.simpleName}_un.gtf
gffread -x -M - -g ${fasta.simpleName}_un.fasta ${gff.simpleName}_un.gtf \
  > ${file_prefix}.fasta
  """
}