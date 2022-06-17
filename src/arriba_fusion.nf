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
      --project [path]                Path to the project folder. Results are saved in this folder.
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3
    
    Input:
      --fastq [path]                  Path to fastq folder.
      --bam [path]                    Path to the bam-containing folder.

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
params.bam = ""
params.fastq = ""
if (params.genome) { params.genome = path(params.genome, checkIfExists: true) } else { exit 1, "No genome specified." }
if (params.gtf) { params.gtf = path(params.gtf, checkIfExists: true) } else { exit 1, "No annotation specified." }

/* Params out */
params.fastp_out = "$params.project/fastp/"
params.index_fasta_out = "$params.project/Indexed_genome/"
params.sort_bam_out = "$params.project/Bam_filtered_sorted/"
params.index_bam_out = "$params.project/Bam_filt_sort_indexed/"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Reference genome : ${params.genome}"
log.info "Genome annotation : ${params.gtf}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

if(params.bam != "") {
    Channel
        .fromPath( params.bam )
        .set { bam_files }
}
else {
    Channel
        .fromFilePairs( params.fastq, size = -1 )
        .set(fastq_files)    
}

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

  //###################### ARRIBA FUSION ########################

  arriba(fastq_files, gtf, genome)

  //################ GRAPHICAL REPRESENTATIONS ##################

}