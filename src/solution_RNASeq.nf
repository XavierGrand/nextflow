nextflow.enable.dsl=2

include { fastp } from "./nf_modules/fastp/main.nf"
include { fasta_from_bed } from "./nf_modules/bedtools/main.nf"
include { index_fasta; mapping_fastq } from './nf_modules/kallisto/main.nf' addParams(mapping_fastq_out: "quantification/")


params.fastq = "data/fastq/*_{1,2}.fastq"

log.info "fastq files: ${params.fastq}"
log.info "fasta file : ${params.fasta}"
log.info "bed file : ${params.bed}"

channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }

channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }
channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .map { it -> [it.simpleName, it]}
  .set { bed_files }

workflow {
  fastp(fastq_files)
  fasta_from_bed(fasta_files, bed_files)
  index_fasta(fasta_from_bed.out.fasta)
  mapping_fastq(index_fasta.out.index.collect(), fastp.out.fastq)
}
