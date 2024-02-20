/*
small pipeline to build a training dataset from whole genome data

input:
- fasta
- fastq
- chromosome
- start position
- stop position

output:
- sort fasta
- sort fastq

example for paired-end data:
./nextflow src/training_dataset.nf -c src/training_dataset.config --fasta "data/genome.fa" --fastq_paired "data/*_R{1,2}.fastq.gz" --chromosome "X" --start 5305683 --stop 5333928 -resume

example for single-end data:
./nextflow src/training_dataset.nf -c src/training_dataset.config --fasta "data/genome.fa" --fastq_single  "data/*_R1.fastq.gz"  --chromosome "X" --start 5305683 --stop 5333928 -resume

*/

params.fastq_paired = ""
params.fastq_single = ""

log.info "fasta files : ${params.fasta}"
log.info "fastq paired files : ${params.fastq_paired}"
log.info "fastq single files : ${params.fastq_single}"
log.info "chromosome : ${params.chromosome}"
log.info "start position : ${params.start}"
log.info "stop position : ${params.stop}"


Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any index files matching: ${params.fasta}" }
  .set { fasta_file }


process build_synthetic_bed {
  tag "${chromosome}:${start}-${stop}"
  cpus 4

  input:
  val chromosome from params.chromosome
  val start from params.start
  val stop from params.stop

  output:
  file "*.bed" into bed_files

  script:
"""
echo "${chromosome}\t${start}\t${stop}" > synthetic.bed
"""
}

process fasta_from_bed {
  tag "${fasta.baseName}"
  cpus 4
  publishDir "results/training/fasta/", mode: 'copy'

  input:
  file fasta from fasta_file
  file bed from bed_files
  val chromosome from params.chromosome

  output:
  file "*.fasta" into fasta_files_extracted

  script:
"""
bedtools getfasta \
-fi ${fasta} -bed ${bed} -fo s${fasta.baseName}.fasta
"""
}

process index_fasta {
  tag "$fasta.baseName"
  cpus 4
  publishDir "results/training/mapping/index/", mode: 'copy'

  input:
    file fasta from fasta_files_extracted

  output:
    file "*.index*" into index_files
    file "*_report.txt" into indexing_report

  script:
"""
bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.baseName}.index &> ${fasta.baseName}_bowtie2_report.txt

if grep -q "Error" ${fasta.baseName}_bowtie2_report.txt; then
  exit 1
fi
"""
}

if ( params.fastq_paired != "" ) {
  Channel
    .fromFilePairs( params.fastq_paired )
    .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq_paired}" }
    .set { fastq_files_paired }

  process mapping_fastq_paired {
    tag "$pair_id"
    cpus 4

    input:
    set pair_id, file(reads) from fastq_files_paired
    file index from index_files.collect()

    output:
    set pair_id, "*.bam" into bam_files_paired
    file "*_report.txt" into mapping_report

    script:
    index_id = index[0]
    for (index_file in index) {
      if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
          index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
      }
    }
  """
  bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
  -1 ${reads[0]} -2 ${reads[1]} 2> \
  ${pair_id}_bowtie2_report.txt | \
  samtools view -Sb - > ${pair_id}.bam

  if grep -q "Error" ${pair_id}_bowtie2_report.txt; then
    exit 1
  fi
  """
  }

  bam_files_paired.into{ bam_files_paired_fa; bam_files_paired_ba}

  process bam_2_fastq_paired {
    tag "$file_id"
    publishDir "results/training/fastq/", mode: 'copy'

    input:
      set file_id, file(bam) from bam_files_paired_fa

    output:
      set file_id, "*.fastq" into fastq_files_extracted
    script:
  """
  samtools fastq -1 s${file_id}_R1.fastq -2 s${file_id}_R2.fastq -F 0x4 ${bam}
  """
  }

  process filter_bam_paired {
    tag "$file_id"
    cpus 4

    input:
      set file_id, file(bam) from bam_files_paired_ba
      file bed from bed_files

    output:
      set file_id, "*.bam" into filtered_bam_files_paired
    script:
  """
  samtools view -@ ${task.cpus} -hb ${bam} -F 0x4 > f${file_id}.bam
  """
  }

  process sort_bam_paired {
    tag "$file_id"
    publishDir "results/training/bams/", mode: 'copy'
    cpus 4

    input:
      set file_id, file(bam) from filtered_bam_files_paired

    output:
      set file_id, "*.bam" into sorted_bam_files_paired

    script:
  """
  samtools sort -@ ${task.cpus} -O BAM -o s${file_id}.bam ${bam}
  """
  }

  process index_bam_paired {
    tag "$file_id"
    publishDir "results/training/bams/", mode: 'copy'

    input:
      set file_id, file(bam) from sorted_bam_files_paired

    output:
      set file_id, "*.bam*" into indexed_bam_file_paired

    script:
  """
  samtools index ${bam}
  """
  }
}


if ( params.fastq_single != "" ) {
  Channel
    .fromPath( params.fastq_single )
    .ifEmpty { error "Cannot find any fastq files matching: ${params.fastq_single}" }
    .map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
    .set { fastq_files_single }

  process mapping_fastq_single {
    tag "$file_id"
    cpus 4

    input:
    set file_id, file(reads) from fastq_files_single
    file index from index_files.collect()

    output:
    set file_id, "*.bam" into bam_files_single
    file "*_report.txt" into mapping_report

    script:
    index_id = index[0]
    for (index_file in index) {
      if (index_file =~ /.*\.1\.bt2/ && !(index_file =~ /.*\.rev\.1\.bt2/)) {
          index_id = ( index_file =~ /(.*)\.1\.bt2/)[0][1]
      }
    }
  """
  bowtie2 --very-sensitive -p ${task.cpus} -x ${index_id} \
  -U ${reads} 2> \
  ${file_id}_bowtie2_report.txt | \
  samtools view -Sb - > ${file_id}.bam

  if grep -q "Error" ${file_id}_bowtie2_report.txt; then
    exit 1
  fi
  """
  }

  bam_files_single.into{ bam_files_single_fa; bam_files_single_ba}

  process bam_2_fastq_single {
    tag "$file_id"

    input:
      set file_id, file(bam) from bam_files_single_fa

    output:
      set file_id, "*.fastq" into fastq_files_extracted
    script:
  """
  samtools fastq -0 s${file_id}.fastq -F 0x4 ${bam}
  """
  }

  process filter_bam_single {
    tag "$file_id"
    cpus 4

    input:
      set file_id, file(bam) from bam_files_single_ba
      file bed from bed_files

    output:
      set file_id, "*.bam" into filtered_bam_files_single
    script:
  """
  samtools view -@ ${task.cpus} -hb ${bam} -F 0x4 > f${file_id}.bam
  """
  }

  process sort_bam_single {
    tag "$file_id"
    publishDir "results/training/bams/", mode: 'copy'
    cpus 4

    input:
      set file_id, file(bam) from filtered_bam_files_single

    output:
      set file_id, "*.bam" into sorted_bam_files_single

    script:
  """
  samtools sort -@ ${task.cpus} -O BAM -o s${file_id}.bam ${bam}
  """
  }

  process index_bam_single {
    tag "$file_id"
    publishDir "results/training/bams/", mode: 'copy'

    input:
      set file_id, file(bam) from sorted_bam_files_single

    output:
      set file_id, "*.bam*" into indexed_bam_file_single

    script:
  """
  samtools index ${bam}
  """
  }
}

