#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      arriba_fusion
========================================================================================================================

star_fusion pipeline :
 * Pipeline dedicated to rna fusion from transcriptomic data using STAR-Fusion software.

Maintainer Xavier Grand <xavier.grand@ens-lyon.fr>

 ****************************************************************
                      Help Message Definition
 ****************************************************************
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

      nextflow ./src/arriba_fusion.nf -c ./src/nextflow.config -profile singularity

    Mandatory arguments:
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3
    
    Input:
      --fastq [path]                  Path to fastq files.
      --bam [path]                    Path to the bam files.

    References:                       Can be downloaded with download_references.sh (not implemented in pipeline).
      --genome [path]                 Path to genome reference fasta file.
      --gtf [path]                    Path to genome annotation gtf file.

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
 
params.bam = ""
params.genome = ""
params.gtf = ""
params.fastq = ""

/* Params out */
params.fastp_out = "Arriba_fastp/"
params.index_fasta_out = "Arriba_Indexed_genome/"
params.sort_bam_out = "Arriba_Bam_filtered_sorted/"
params.index_bam_out = "Arriba_Bam_filt_sort_indexed/"
params.arriba_out = "Arriba_Arriba_results/"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Reference genome : ${params.genome}"
log.info "Genome annotation : ${params.gtf}"
if(params.bam != "") {
  bam_list = "${params.bam}/*.bam"
  log.info "bam files (--bam): ${bam_list}"
}
else {
  log.info "fastq files (--fastq): ${params.fastq}"
}

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

if(params.bam != "") {
    Channel
        .fromPath( bam_list )
        .ifEmpty { error "Cannot find any bam files in: ${params.bam}" }
        .map { it -> [it.simpleName, it]}
        .set { bam_files }
}
else {
    Channel
        .fromFilePairs( params.fastq, size: -1)
        .set { fastq_files }
}

Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any fasta files in: ${params.genome}" }
  .map { it -> [it.simpleName, it]}
  .set { genome }

Channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any annotation files in: ${params.gtf}" }
  .map { it -> [it.simpleName, it]}
  .set { gtf }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { fastp } from './nf_modules/fastp/main.nf'
include { fastqc_fastq as fastqc_raw } from './nf_modules/fastqc/main.nf' addParams(fastqc_fastq_out: "$params.project/01_fastqc_raw/")
include { fastqc_fastq as fastqc_preprocessed } from './nf_modules/fastqc/main.nf' addParams(fastqc_fastq_out: "$params.project/02_fastqc_preprocessed/")
include { multiqc } from './nf_modules/multiqc/main.nf' addParams(multiqc_out: "$params.project/QC/")
include { index_with_gtf } from './nf_modules/star/main_2.7.8a.nf' addParams(star_mapping_fastq_out: "$params.project/STAR_index/")
include { mapping_fastq_withChimeric } from './nf_modules/star/main_2.7.8a.nf' addParams(star_mapping_fastq_out: "$params.project/STAR/")
include { arriba } from "./nf_modules/arriba/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  if(params.bam == ""){
    fastp(fastq_files)
    // fastqc_raw(fastq_files.collect())
    // fastqc_preprocessed(fastp_out.fastq.collect())
    // multiqc(fastqc_raw_out.report)
    // .mix(
    //   fastqc_preprocessed.out.report
    //   ).collect()
    index_with_gtf(genome, gtf)
    // mapping_fastq_withChimeric(index_fasta_out.index, fastp_out.fastq)
    // filter_bam_quality(mapping_fastq_withChimeric_out.bam)
    // arriba()
  }
  else {
    arriba(bam_files, gtf, genome)
  }

}