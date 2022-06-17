#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      star-fusion
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

      nextflow ./src/star_fusion.nf -c ./src/nextflow.config -profile singularity

    Mandatory arguments:
      --project [path]                Path to the project folder containing fastq folder. Results are saved in this folder.
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3

    References:
      --lib [path]                    Path to CTAT resource library.
                                      Available at https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/

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
params.fastq = "${params.project}/fastq"

/*
if (params.genome)          { params.genome = path(params.genome, checkIfExists: true) } else { exit 1, "No genome specified." }
*/

/* Params out */
params.star_fusion_out = "$params.project/predicted_fusion/"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Genome library folder : ${params.lib}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

Channel
  .fromFilePairs( params.fastq, size: -1 )
  .set { fastq_files }

Channel
  .fromPath( params.lib )
  .set { lib_folder }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { star_fusion } from "./nf_modules/star-fusion/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  //######################## STAR FUSION ########################

  star_fusion(lib_folder, fastq_files)

  //################ GRAPHICAL REPRESENTATIONS ##################

}