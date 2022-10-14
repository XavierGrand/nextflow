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

      nextflow ./src/Readthrough_analysis.nf -c ./src/nextflow.config -profile singularity

    Configuration argument:
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3

    Samples:
      -fastq [path]                   Path to fastq folder.

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
 
params.fastp_out = "fastp/"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Annotation: ${params.gtf}"
log.info "Genome fasta file: ${params.fasta}"
log.info "Genome index location: ${params.idx}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

/*
Channel
  .fromFilePairs( params.fastq, size: -1 )
  .set { fastq_files }

Channel
  .fromPath( params.gtf )
  .set { gtf_file }
*/

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

fastqc_mod = "./nf_modules/fastqc/main.nf"
include { fastqc_fastq as fastqc_raw } from fastqc_mod addParams(fastqc_fastq_out: "01_fastqc_raw/")
include { fastqc_fastq as fastqc_preprocessed } from fastqc_mod addParams(fastqc_fastq_out: "02_fastqc_preprocessed/")
include { multiqc } from './nf_modules/multiqc/main.nf' addParams(multiqc_out: "QC/")
include { fastp } from "./nf_modules/fastp/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

}
