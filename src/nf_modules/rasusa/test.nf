nextflow.enable.dsl=2

/*
./nextflow src/nf_modules/rasusa/test.nf -c src/nextflow.config -profile docker --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --fastq "data/tiny_dataset/fastq/tiny_R1.fastq"
./nextflow src/nf_modules/rasusa/test.nf -c src/nextflow.config -profile docker --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" --coverage 1.0
./nextflow src/nf_modules/rasusa/test.nf -c src/nextflow.config -profile docker --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --fastq "data/tiny_dataset/fastq/tiny_R1.fastq" --size "1Mb"
*/

params.fastq = "data/fastq/*R{1,2}*"
params.fasta = "data/fasta/*.fasta"
params.coverage = ""
params.size = ""

include { sample_fastq } from "./main.nf" addParams(sample_fastq_coverage: params.coverage, sample_fastq_size: params.size, sample_fastq_out: "sample/")

channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }

channel
  .fromPath( params.fasta )
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }

workflow {
  sample_fastq(fastq_files, fasta_files.collect())
}