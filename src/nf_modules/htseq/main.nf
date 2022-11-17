
version = "0.13.5"
container_url = "lbmc/htseq:${version}"

params.htseq_out = ""



process gff3_2_gtf {
  container = "dceoy/cufflinks"
  label "small_mem_mono_cpus"

    input:
        tuple val(genome_id), path(gff3_file)
    output:
        path "${genome_id}.gtf", emit: gtf
    script:
"""
gffread ${gff3_file} -T -o ${genome_id}.gtf
"""
}


params.htseq_param = "yes"
process htseq_count {
    container = "${container_url}"
    label "big_mem_mono_cpus"
    tag "file_id: $file_id"
    if (params.htseq_out != "") {
        publishDir "results/${params.htseq_out}", mode: 'copy'
    }
    input:
      tuple val(file_id), path(bam), path(bai)
      path(gtf)

    output:
      path "${file_id}.tsv", emit: counts

  script:
"""
htseq-count -r pos -a 10 -s ${params.htseq_param} -t exon -i gene_id $bam $gtf > ${file_id}.tsv
"""
}

workflow htseq_count_with_gff {
  take:
    bam_tuple
    gff_file
  main:
    gff3_2_gtf(gff_file)
    htseq_count(bam_tuple,gff3_2_gtf.out.gtf)
  emit:
    counts = htseq_count.out.counts
}
