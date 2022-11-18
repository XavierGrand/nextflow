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
      --fastq [path]                  Path to fastq folder.
      --bam [path]                    Path to the bam folder (indexed and sorted).

    References:                       Can be downloaded with download_references.sh (not implemented in pipeline).
      --genome [path]                 Path to genome reference fasta file.
      --index [path]                  Path to STAR-indexed genome folder.
      --gtf [path]                    Path to genome annotation gtf file.

    Optionnal:
      --arriba_options [str]          Tunable options for arriba, among: [-c Chimeric.out.sam] [-b blacklists.tsv] 
                                      [-k known_fusions.tsv] [-d structural_variants_from_WGS.tsv] [-t tags.tsv] 
                                      [-p protein_domains.gff3] [OPTIONS]
                                      Options have to be given between quotes, i.e.: 
                                      --arriba_options "-c Chimeric.out.sam -b blacklists.tsv"
                                      Except: -x Aligned.out.sam -g annotation.gtf -a assembly.fa -o fusions.tsv 
                                      -O fusions.discarded.tsv

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
params.fastp_out = "02_fastp"
params.index_fasta_out = "04_Indexed_genome"
params.star_mapping2fusion_fastq_out = "06_mapping2fusion"
params.sort_bam_out = "07_sort_bam"
params.index_bam_out = "08_index_bam"
params.arriba_out = "09_Arriba_results"

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Reference genome : ${params.genome}"
log.info "Genome annotation : ${params.gtf}"
if(params.bam != "") {
  bam_list = "${params.bam}/*_aligned_sorted.bam"
  log.info "bam files (--bam): ${bam_list}"
}
else {
  fastq_list = "${params.fastq}/*_{1,2}.fastq"
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
include { fastqc_fastq as fastqc_raw } from './nf_modules/fastqc/main.nf'
include { fastqc_fastq as fastqc_preprocessed } from './nf_modules/fastqc/main.nf'
include { multiqc } from './nf_modules/multiqc/main.nf'
include { index_with_gtf } from './nf_modules/star/main_2.7.8a.nf'
include { mapping_fastq_withChimeric } from './nf_modules/star/main_2.7.8a.nf'
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