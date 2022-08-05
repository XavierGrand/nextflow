version = "2.7.8a"
container_url = "xgrand/star:${version}"

params.star_mapping_fastq_out = ""


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


process index_with_gtf {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  if (params.star_index_out != "") {
    publishDir "results/${params.star_index_out}", mode: 'copy'
  }

  input:
    tuple val(genome_id), path(genome_fasta)
    path gtf_file

  output:
    tuple val(genome_id), path ("*"), emit: index

  script:
"""
STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles ${genome_fasta}  \
--sjdbGTFfile ${gtf_file} \
--genomeSAindexNbases 13 # min(14, log2(GenomeLength)/2 - 1)
"""
}

workflow index_with_gff {
  take:
    genome_fasta
    gff_file
  main:
    gff3_2_gtf(gff_file)
    index_with_gtf(genome_fasta,gff3_2_gtf.out.gtf)
  emit:
    report = index_with_gtf.out.index
}


process index_without_gff {
  container = "${container_url}"
  label "big_mem_multi_cpus"

  input:
    tuple val(genome_id), path(genome_fasta)

  output:
    tuple val(genome_id), path ("*"), emit: index

  script:
"""
STAR --runThreadN ${task.cpus} --runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles ${genome_fasta}  \
--genomeSAindexNbases 13 # min(14, log2(GenomeLength)/2 - 1)
"""
}


process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  if (params.star_mapping_fastq_out != "") {
    publishDir "results/${params.star_mapping_fastq_out}", mode: 'copy'
  }

  input:
    tuple val(index_id), path(index)
    tuple val(reads_id), path(reads) 

  output:
    path "*.Log.final.out", emit: report
    tuple val(reads_id), path("*.bam"), emit: bam

  script:
if (reads_id instanceof List){
    file_prefix = reads_id[0]
  } else {
    file_prefix = reads_id
  }

if (reads.size() == 2)
"""
mkdir -p index
mv ${index} index/
STAR --runThreadN ${task.cpus} \
--genomeDir index/ \
--readFilesCommand zcat \
--readFilesIn ${reads[0]} ${reads[1]} \
--outFileNamePrefix ${reads_id}. \
--alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif

mv ${reads_id}.Aligned.sortedByCoord.out.bam ${reads_id}.bam
"""
else
"""
mkdir -p index
mv ${index} index/
STAR --runThreadN ${task.cpus} \
--genomeDir index/ \
--readFilesCommand zcat \
--readFilesIn ${reads} \
--outFileNamePrefix ${reads_id}. \
--alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif

mv ${reads_id}.Aligned.sortedByCoord.out.bam ${reads_id}.bam
"""
}

process mapping_fastq_withChimeric {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  if (params.star_mapping_fastq_out != "") {
    publishDir "results/${params.star_mapping_fastq_out}", mode: 'copy'
  }

  input:
    tuple val(index_id), path(index)
    tuple val(reads_id), path(reads) 

  output:
    path "*.Log.final.out", emit: report
    tuple val(reads_id), path("*.bam"), emit: bam

  script:
if (reads_id instanceof List){
    file_prefix = reads_id[0]
  } else {
    file_prefix = reads_id
  }

if (reads.size() == 2)
"""
mkdir -p index
mv ${index} index/
STAR --runThreadN ${task.cpus} \
--genomeDir index/ \
--readFilesCommand zcat \
--readFilesIn ${reads[0]} ${reads[1]} \
--outFileNamePrefix ${reads_id}. \
--alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--chimOutType WithinBAM

mv ${reads_id}.Aligned.sortedByCoord.out.bam ${reads_id}.bam
"""
else
"""
mkdir -p index
mv ${index} index/
STAR --runThreadN ${task.cpus} \
--genomeDir index/ \
--readFilesCommand zcat \
--readFilesIn ${reads} \
--outFileNamePrefix ${reads_id}. \
--alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--chimOutType WithinBAM

mv ${reads_id}.Aligned.sortedByCoord.out.bam ${reads_id}.bam
"""
}

process mapping_withindex {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  if (params.star_mapping_fastq_out != "") {
    publishDir "results/${params.star_mapping_fastq_out}", mode: 'copy'
  }

  input:
    path(index)
    tuple val(reads_id), path(reads) 

  output:
    path "*.Log.final.out", emit: report
    tuple val(reads_id), path("*.bam"), emit: bam

  script:
if (reads_id instanceof List){
    file_prefix = reads_id[0]
  } else {
    file_prefix = reads_id
  }

if (reads.size() == 2)
"""
STAR --runThreadN ${task.cpus} \
--genomeDir ${index} \
--readFilesCommand zcat \
--readFilesIn ${reads[0]} ${reads[1]} \
--outFileNamePrefix ${reads_id}. \
--alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif

mv ${reads_id}.Aligned.sortedByCoord.out.bam ${reads_id}.bam
"""
else
"""
STAR --runThreadN ${task.cpus} \
--genomeDir ${index} \
--readFilesCommand zcat \
--readFilesIn ${reads} \
--outFileNamePrefix ${reads_id}. \
--alignIntronMax 10000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif

mv ${reads_id}.Aligned.sortedByCoord.out.bam ${reads_id}.bam
"""
}