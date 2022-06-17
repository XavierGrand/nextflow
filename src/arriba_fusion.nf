#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      arriba_fusion
========================================================================================================================

star_fusion pipeline :
 * Pipeline dedicated to rna fusion from transcriptomic data using STAR-Fusion software.

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
      --project [path]                Path to the project folder containing fastq folder. Results are saved in this folder.
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3

    References:
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
 
params.project = ""
params.paired = true
params.fastq = "${params.project}/fastq/*_R{1,2}.fastq*"

/*
if (params.genome)          { params.genome = path(params.genome, checkIfExists: true) } else { exit 1, "No genome specified." }
*/

/* Params out */

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Reference genome : ${params.genome}"
log.info "Genome annotation : ${params.genome}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

Channel
  .fromFilePairs( params.fastq, size: -1 )
  .set { fastq_files }

Channel
  .fromPath( params.genome )
  .set { genome }

Channel
  .fromPath( params.gtf )
  .set { gtf }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { arriba } from "./nf_modules/arriba/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  fastq_files.view()

  //###################### ARRIBA FUSION ########################

  arriba(fastq_files, gtf, genome)

  //################ GRAPHICAL REPRESENTATIONS ##################

}