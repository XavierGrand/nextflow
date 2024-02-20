version = "1.0"
container_url = "lbmc/bioawk:${version}"

params.fasta_to_transcripts_lengths = ""
params.fasta_to_transcripts_lengths_out = ""
process fasta_to_transcripts_lengths {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.fasta_to_transcripts_lengths_out != "") {
    publishDir "results/${params.fasta_to_transcripts_lengths_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("${fasta.simpleName}_transcripts_lengths.tsv"), emit: tsv

  script:
"""
bioawk -c fastx '{print(\$name" "length(\$seq))}' ${fasta} > ${fasta.simpleName}_transcripts_lengths.tsv
"""
}