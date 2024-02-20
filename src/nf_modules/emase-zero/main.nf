version = "0.3.1"
container_url = "lbmc/emase-zero:${version}"

include { g2tr } from "./../kb/main.nf"
include { bam2ec } from "./../alntools/main.nf"
include { fasta_to_transcripts_lengths } from "./../bioawk/main.nf"


params.count = "-m 2"
params.count_out = ""
workflow count {
  take:
    bam_idx
    fasta
    gtf

  main:
    g2tr(gtf)
    fasta_to_transcripts_lengths(fasta)
    bam2ec(bam_idx, fasta_to_transcripts_lengths.out.tsv.collect())
    emase(bam2ec.out.bin, fasta.collect(), bam2ec.out.tsv, g2tr.out.g2t.collect())

  emit:
    count = emase.out.count
}

process emase {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.count_out != "") {
    publishDir "results/${params.count_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bin)
    tuple val(fasta_id), path(fasta)
    tuple val(transcript_length_id), path(transcript_length)
    tuple val(gene_to_transcript_id), path(gene_to_transcript)

  output:
    tuple val(file_id), path("${bin.simpleName}.quantified*"), emit: count
    path "*_report.txt", emit: report

  script:
"""
grep ">" ${fasta} | sed 's/>//' > tr_list.txt
emase-zero ${params.count} \
  -o ${bin.simpleName}.quantified \
  -l ${transcript_length} \
  -g ${gene_to_transcript} \
  ${bin} &> ${file_id}_emase-zero_report.txt
"""
}