version = "4.6.1.0"
container_url = "xgrand/gatk:${version}"


params.pathseq_out = ""
process pathseq {
    container = "${container_url}"
  label "huge_mem_multi_cpus"
  tag "$file_id"
  if (params.pathseq_out != "") {
    publishDir "results/${params.pathseq_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam), path(bai)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf

  script:
  memory = "${task.memory}" - ~/\s*GB/
"""
gatk --java-options "-Xmx${memory}g" PathSeqPipelineSpark \
	    --input ${file} \
	    --filter-bwa-image ${pathseqdb}/pathseq_host.fa.img \
	    --kmer-file ${pathseqdb}/pathseq_host.bfi \
	    --min-clipped-read-length 31 \
	    --microbe-bwa-image ${pathseqdb}/pathseq_microbe.fa.img \
	    --microbe-fasta ${pathseqdb}/pathseq_microbe.fa \
	    --taxonomy-file ${pathseqdb}/pathseq_taxonomy.db \
	    --output ${outpath}/${samplename}.pathseq.complete.bam \
	    --scores-output ${outpath}/${samplename}.pathseq.complete.csv \
	    --is-host-aligned false \
	    --filter-duplicates false \
	    --min-score-identity .7
"""
}