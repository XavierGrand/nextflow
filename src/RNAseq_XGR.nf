#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      RNAseq_XGR
========================================================================================================================

RNAseq_XGR pipeline :
 * Pipeline dedicated to transcriptomic analysis.

Maintainer Xavier Grand <xavier.grand@ens-lyon.fr>

 ****************************************************************
                      Help Message Definition
 ****************************************************************
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

      nextflow ./src/RNAseq_XGR.nf -c ./src/nextflow.config -profile singularity

    Mandatory arguments:
      --project [path]                Path to the project folder. Results are saved in this folder.
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3

    References:
      --fasta [path]                  Path to genome fasta file.
      --gtf [path]                    Path to the gtf annotation file.
      --idx [path]                    Path to the STAR indexed genome (optional). (If allready computed)

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
 
/* Arguments */
project = params.project
params.fastq = "${project}/fastq/*_{1,2}.fq.gz"
params.gtf = ""
params.fasta = ""
params.idx = ""

params.fastp_out = "$params.project/fastp/"
params.star_mapping_fastq_out = "$params.project/STAR/"
params.star_index_out = "$params.project/STARindex/"

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

Channel
  .fromFilePairs( params.fastq, size: -1 )
  .set { fastq_files }

Channel
  .fromPath( params.gtf )
  .set { gtf_file }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

fastqc_mod = "./nf_modules/fastqc/main.nf"
include { fastqc_fastq as fastqc_raw } from fastqc_mod addParams(fastqc_fastq_out: "$params.project/01_fastqc_raw/")
include { fastqc_fastq as fastqc_preprocessed } from fastqc_mod addParams(fastqc_fastq_out: "$params.project/02_fastqc_preprocessed/")
include { multiqc } from './nf_modules/multiqc/main.nf' addParams(multiqc_out: "$params.project/QC/")
include { fastp } from "./nf_modules/fastp/main.nf"
include { index_with_gtf } from "./nf_modules/star/main_2.7.8a.nf"
include { mapping_fastq } from "./nf_modules/star/main_2.7.8a.nf"
include { mapping_withindex } from "./nf_modules/star/main_2.7.8a.nf"
include { htseq_count } from "./nf_modules/htseq/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  //########################## PREPROCESSING ####################   
  // fastp
  fastp(fastq_files)

  //########################## QUALITY CHECKS ###################

  // fastqc_rawdata
  fastqc_raw(fastq_files)
  // fastqc_processed
  fastqc_preprocessed(fastp.out.fastq.map { it -> [it [0], it[1]]})
  // multiqc
  multiqc(
    fastqc_raw.out.report
    .mix(
      fastqc_preprocessed.out.report
      ).collect()
  )

  //############ GENOME INDEXATION AND MAPPING ###################

  if (params.idx == "") {
    Channel
      .fromPath( params.fasta )
      .ifEmpty { error "Cannot find any files matching: ${params.fasta}" }
      .map{it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
      .set { genome_file }
    
    index_with_gtf(genome_file, gtf_file.collect())
    mapping_fastq(index_with_gtf.out.index.collect(), fastp.out.fastq)
    htseq_count(mapping_fastq.out.bam, gtf_file)
  }
  else {
    idx_genome = "${params.idx}"
    Channel
      .of( idx_genome )
      .set { genome_indexed_input }
    genome_indexed_input.view()
    mapping_withindex(fastp.out.fastq)
    htseq_count(mapping_withindex.out.bam, gtf_file)
  }
}
