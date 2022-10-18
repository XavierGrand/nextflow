#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      Readthrough analysis
========================================================================================================================

Readthrough analysis pipeline :
 * Pipeline dedicated to readthrough analysis.

Maintainer Xavier Grand <xavier.grand@ens-lyon.fr>

 ****************************************************************
                      Help Message Definition
 ****************************************************************
*/

def helpMessage() {
    log.info"""
    Usage:
    Pipeline dedicated to readthrough analysis of short-reads paired-ends RNAseq.
    The typical command for running the pipeline is as follows:

      nextflow ./src/readthrough_analysis.nf -c ./src/nextflow.config -profile singularity

    Configuration argument:
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3

    Samples:
      --fastq [path]                  Path to fastq folder.
      --bam [path]                    Path to bam folder.

    References:
      --fasta [path]                  Path to genome fasta file.
      --gtf [path]                    Path to the gtf annotation file.

    Readthrough:
      --rt_length [int]               Length of the readthrough range to consider (default=10000).

    Help:                             Display this help message.
      --help
      --h
    
    """.stripIndent()
}

// Show help message

params.help = ""
params.h = ""

if (params.help || params.h) {
    helpMessage()
    exit 0
}

/*
 ****************************************************************
                      Default Parameters
 ****************************************************************
*/

// params.fastq = "*{1,2}.fastq"
params.bam = "data/*.bam"
params.gtf = ""
params.fasta = ""
params.flag_fwd = "99,147"
params.flag_rev = "83,163"

params.rtranger_out = "05_rtranger/"
params.split_rt_bam_out = "06_splited_bam/"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

// log.info "Genome gtf file: ${params.gtf}"
// log.info "Genome fasta file: ${params.fasta}"
log.info "Bam files: ${params.bam}"
log.info "FLAG value for Forward read pairs: ${params.flag_fwd}"
log.info "FLAG value for Reverse read pairs: ${params.flag_rev}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

/*
if (params.fastq != "") {
  Channel
  .fromFilePairs( params.fastq, size: -1 )
  .set { fastq_files }
}
*/

Channel
  .fromPath( params.bam )
  .ifEmpty { error "Cannot find any genome files matching: ${params.bam}" }
  .map( it -> [it.baseName, it])
  .set { bam_input }

Channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any gtf files matching: ${params.gtf}" }
  .set { gtf_file }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

/*
fastqc_mod = "./nf_modules/fastqc/main.nf"
include { fastqc_fastq as fastqc_raw } from fastqc_mod addParams(fastqc_fastq_out: "01_fastqc_raw/")
include { fastqc_fastq as fastqc_preprocessed } from fastqc_mod addParams(fastqc_fastq_out: "02_fastqc_preprocessed/")
include { multiqc } from './nf_modules/multiqc/main.nf' addParams(multiqc_out: "03_MultiQC/")
include { fastp } from "./nf_modules/fastp/main.nf" addParams(fastp_out = "04_fastp/")
*/
include { rtranger } from "./nf_modules/rtranger/main.nf"
include { split_rt_bam } from "./nf_modules/samtools/main.nf"
include { index_bam as index_bam_fwd } from "./nf_modules/samtools/main.nf" addParams(index_bam_out: "06_splited_bam/")
include { index_bam as index_bam_rev } from "./nf_modules/samtools/main.nf" addParams(index_bam_out: "06_splited_bam/")

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {
  rtranger(gtf_file)
  split_rt_bam(bam_input)
  index_bam_fwd(split_rt_bam.out.bam_forward)
  index_bam_rev(split_rt_bam.out.bam_reverse)

}
