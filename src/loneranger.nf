#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                              Sataus: In Development...
========================================================================================================================
*/

/*
========================================================================================================================
                                                      LoneRanger
========================================================================================================================

LoneRanger pipeline :
* Pipeline dedicated to detect microbes in scRNAseq 10X with Cellranger, from Galeano Niño et al. 2023.

Maintainer Xavier Grand <xavier.grand@inserm.fr>

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
      --input [path]                  Path to fastq folder.
      --design [tsv]                  Path to design tsv file.

    References:
      --genome [path]                 Path to genome reference fasta file.
      --index [path]                  Path to STAR-indexed genome folder.
      --gtf [path]                    Path to genome annotation gtf file.

    Optionnal:

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
 
params.sample_list = ["HBV","HDV","Mock","HBV-HDV"]
params.input = "./data/fastq"
params.fastq = "${params.input}/*_{R1,R2,I1}*.fastq.gz"
params.genome = "/home/xavier/Data/Genome/10X_HumanReference/refdata-gex-GRCh38-2024-A"

/* Params out */
params.cellranger_count_out = "01_cellranger_output"
params.pathseq_out = "02_pathseq"
params.invadeseq_matrix_out = "03_matrix"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Reference genome : ${params.genome}"
log.info "Input directory : ${params.input}"
log.info "Loaded fastq files (--fastq): ${params.fastq}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

Channel
  .fromList ( params.sample_list )
  .map { v -> v }
  .flatten()
  .set { sample }

Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any genome files matching: ${params.genome}" }
  .set { refdata }

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { cellranger_count } from './nf_modules/cellranger/9.0.0/main.nf'
include { pathseq } from "./nf_modules/gatk/4.6.1.0/main.nf"
//include { invadeseq_matrix } from "./nf_modules/invadeseq/1.0/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

/*
 ****************************************************************
                      Cellranger
 ****************************************************************
*/
  
  cellranger_count(sample, params.input, refdata)

/*
 ****************************************************************
                          Pathseq
 ****************************************************************
*/



/*
 ****************************************************************
                    invadeseq
 ****************************************************************
*/

}