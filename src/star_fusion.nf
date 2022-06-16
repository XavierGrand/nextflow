#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      star_fusion
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
      --project [path]                Path to the project folder. Results are saved in this folder.
      --input [path]                  Path to the folder containing fast5 files. If skip basecalling option enabled, path to fastq files folder.
      -profile [str]                  Configuration profile to use.
                                      Available: docker, singularity, podman, psmn, ccin2p3

    References:                       If not specified, use of default transcriptome <EBI ID>.
      --genome [file]                 Path to the fasta file containing the genome.
      --gtf [file]                    Path to the gtf file containing the genome annotation.

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

/*
if (params.transcriptome)   { params.transcriptome = path(params.transcriptome, checkIfExists: true) } else { exit 1, "No transcriptome specified." }
if (params.genome)          { params.genome = path(params.genome, checkIfExists: true) } else { exit 1, "No genome specified." }
if (params.gtf)             { params.gtf = path(params.gtf, checkIfExists: true) } else { exit 1, "No annotation specified." }
if (params.input)           { params.input = path(params.input, checkIfExists: true) } else { exit 1, "$params.input does not exists." }
if (params.mode != "gpu")   { params.mode = "cpu" } 
*/

/* Params out */

params.basecalling_out = "$params.project/Basecalling/"
params.fastq_out = "$params.project/fastq/"
params.ref_out = "$params.project/Ref/"
params.pycoQC_out = "$params.project/pycoQC/"
params.porechop_out = "$params.project/porechop/"
params.cutadapt_out = "$params.project/cutadapt/"
params.minimap2_transcriptome_out = "$params.project/minimap2_transcriptome/"
params.minimap2_genome_out = "$params.project/minimap2_genome/"
params.sort_bam_out = "$params.project/minimap2_genome/"
params.salmon_out = "$params.project/salmon/"
params.stringtie_out = "$params.project/stringtie2/"
params.start_position_counts_out = "$params.project/start_positions/"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "fast5/q folder : ${params.input}"
log.info "5'RACE adapter sequence : ${params.adapt}"
log.info "Guppy basecalling calculation mode : ${params.mode}"
log.info "Transcriptome file : ${params.transcriptome}"
log.info "Genome file : ${params.genome}"
log.info "Genome annotation file : ${params.gtf}"

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

Channel
    .of( params.input )
    .ifEmpty { error "No fast5/q folder defined." }
    .set { input }

Channel
    .fromPath( params.transcriptome )
    .ifEmpty { error "No transcriptome defined, a fasta file containing transcripts linked to 5pRACE adapter." }
    .set { transcriptome }
