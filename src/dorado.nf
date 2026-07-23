#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                Dorado
========================================================================================================================

Banjul analysis pipeline :
 * Pipeline dedicated to Nanopore basecalling.

Author Xavier Grand <xavier.grand@inserm.fr>

Creation date: 23/07/2026 by Xavier Grand

 ****************************************************************
                      Help Message Definition
 ****************************************************************
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

      nextflow ./src/bolero.nf -c ./src/nextflow.config -profile singularity --input <PATH_TO_INPUT_FOLDER>

    Nextflow parameters:
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3, etc.
                                      User can use his own, associated with a nextflow configuration file.
                                      Refer to https://www.nextflow.io/docs/latest/config.html#configuration-file.

    Inputs arguments:
      --input [path]                  Path to the folder containing fast5 files. 
                                      If skip basecalling option enabled, path to fastq files folder.

    Nanopore basecalling:
      --kit_barcoding [str]           Nanopore barcoding kit. Default = "EXP-PBC001"
      --model [str]                   Model to basecalling with dorado. Available: fast, hac, sup.

    Optional basecalling parameters:
      --min_qscore [int]              Minimum quality score threshold, default = 7.0.

    Help:                           Display this help message.
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

params.input = ""
params.kit_barcoding = "EXP-PBC001"
params.model = "sup"
params.min_qscore = 7.0

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Input path: ${params.input}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

Channel
    .of( params.input )
    .ifEmpty { error "No fast5/q/pod5 folder defined." }
    .set { input }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { pod5convert } from "./nf_modules/dorado/1.1.1/main.nf"
include { basecall_bc } from "./nf_modules/dorado/1.1.1/main.nf"
include { demux } from "./nf_modules/dorado/1.1.1/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

    //####################### LOADING DATA ########################
    //####################### & BASECALLING #######################

    //Multiplexed samples:
    if ( params.kit_barcoding != "" ) {
        //pod5convert(input)
        //basecall_bc(pod5convert.out.pod5)
        // Modif avant implementation input as pod5:
        basecall_bc(input)
        demux(basecall_bc.out.basecalled)
        demux.out.fastq.flatten()
            .map{it -> [it.getSimpleName(), it]}
            .set{tuples_fastq}
    }
    //Single sample:
    else if ( params.kit_barcoding == "" ) {
        pod5convert(input)
        params.dorado_parameters = "--emit-fastq"
        basecall(pod5convert.out.pod5)
        basecall.out.basecalled.map{it -> [it.getSimpleName(), it]}
            .set{tuples_fastq}
    }
}
