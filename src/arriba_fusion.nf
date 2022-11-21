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
 
params.fastq = ""
params.bam = ""
params.genome = ""
params.gtf = ""

/* Params out */
params.fastp_out = "02_fastp"
params.star_index_out = "04_Indexed_genome"
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
if(params.bam) {
  bam_list = "${params.bam}/*_aligned_sorted.bam"
  log.info "Loaded bam files (--bam): ${bam_list}"
}
else {
  fastq_list = "${params.fastq}/*_{R1,R2}.fastq.gz"
  log.info "Loaded fastq files (--fastq): ${fastq_list}"
}

/*
 ****************************************************************
                        Channel definitions
 ****************************************************************
*/

Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any genome files matching: ${params.genome}" }
  .map( it -> [it.baseName, it])
  .set { genome_file }

Channel
    .fromPath( params.gtf )
    .ifEmpty { error "Cannot find any gtf files matching: ${params.gtf}" }
    .map( it -> [it.baseName, it])
    .set { gtf_file }

if(params.bam) {
    Channel
        .fromPath( bam_list )
        .ifEmpty { error "Cannot find any bam files matching: ${params.bam_file}" }
        .map { it -> [it.simpleName, it]}
        .set { bam_files }
}
else {
    Channel
        .fromFilePairs( fastq_list, size: -1)
        .set { fastq_files }
}

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
include { mapping2fusion } from './nf_modules/star/main_2.7.8a.nf'
include { arriba } from "./nf_modules/arriba/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  if(params.fastq != ""){
    fastp(fastq_files)
    fastqc_raw(fastq_files.collect())
    fastqc_preprocessed(fastp.out.fastq.collect())
    /* multiqc(fastqc_raw.out.report)
     .mix(
       fastqc_preprocessed.out.report
       ).collect()
    */
    index_with_gtf(genome_file, gtf_file)
    mapping2fusion(index_with_gtf.out.index, fastp.out.fastq)
    filter_bam_quality(mapping2fusion.out.bam)
    arriba(filter_bam_quality.out.bam, gtf_file, genome_file)
  }
  else {
    arriba(bam_files, gtf_file, genome_file)
  }
}