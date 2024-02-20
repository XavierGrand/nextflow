#! /bin/sh
./nextflow src/nf_modules/rasusa/test.nf -c src/nextflow.config -profile docker --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --fastq "data/tiny_dataset/fastq/tiny_R1.fastq"
./nextflow src/nf_modules/rasusa/test.nf -c src/nextflow.config -profile docker --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --fastq "data/tiny_dataset/fastq/tiny_R{1,2}.fastq" --coverage 1.0
./nextflow src/nf_modules/rasusa/test.nf -c src/nextflow.config -profile docker --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --fastq "data/tiny_dataset/fastq/tiny_R1.fastq" --size "1Mb"