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
    Pipeline dedicated to transcriptomic analysis of short-reads paired-ends RNAseq.
    The typical command for running the pipeline is as follows:

      nextflow ./src/RNAseq_XGR.nf -c ./src/nextflow.config -profile singularity

    Mandatory arguments:
      --project [path]                Path to the fastq folder.
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
params.project = ""
params.fastq = "${params.project}/fastq/*{1,2}.fq.gz"
params.gtf = ""
params.fasta = ""
params.idx = ""
params.filter_bam_mapped = "-F 268 -f 1 -q 10"

params.fastp_out = "fastp/"
params.star_mapping_fastq_out = "STAR/"
params.star_index_out = "STARindex/"
params.sort_bam_out = "STAR/"
params.filter_bam_out = "STAR/"
params.htseq_out = "HTseq/"

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
include { fastqc_fastq as fastqc_raw } from fastqc_mod addParams(fastqc_fastq_out: "01_fastqc_raw/")
include { fastqc_fastq as fastqc_preprocessed } from fastqc_mod addParams(fastqc_fastq_out: "02_fastqc_preprocessed/")
include { multiqc } from './nf_modules/multiqc/main.nf' addParams(multiqc_out: "QC/")
include { fastp } from "./nf_modules/fastp/main.nf"
include { index_with_gtf } from "./nf_modules/star/main_2.7.8a.nf"
include { mapping_fastq } from "./nf_modules/star/main_2.7.8a.nf"
include { mapping_withindex } from "./nf_modules/star/main_2.7.8a.nf"
include { htseq_count } from "./nf_modules/htseq/main.nf"

include { filter_bam_mapped } from "./nf_modules/samtools/main.nf"
include { stats_bam } from "./nf_modules/samtools/main.nf"
include { sort_bam } from "./nf_modules/samtools/main.nf"
include { index_bam } from "./nf_modules/samtools/main.nf"

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  //########################## PREPROCESSING ####################   
  // fastp
  fastp(fastq_files)

  //############ GENOME INDEXATION AND MAPPING ###################

  if (params.idx == "") {
    Channel
      .fromPath( params.fasta )
      .ifEmpty { error "Cannot find any files matching: ${params.fasta}" }
      .map{it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
      .set { genome_file }
    
    index_with_gtf(genome_file, gtf_file.collect())
    mapping_fastq(index_with_gtf.out.index.collect(), fastp.out.fastq)
    filter_bam_mapped(mapping_fastq.out.bam)
  }
  else {
    idx_genome = "${params.idx}"
    ref_1  = Channel.fromPath(params.idx +'/chrStart.txt')
	  ref_2  = Channel.fromPath(params.idx +'/chrNameLength.txt')
	  ref_3  = Channel.fromPath(params.idx +'/chrName.txt')
	  ref_4  = Channel.fromPath(params.idx +'/chrLength.txt')
	  ref_5  = Channel.fromPath(params.idx +'/exonGeTrInfo.tab')
	  ref_6  = Channel.fromPath(params.idx +'/exonInfo.tab')
	  ref_7  = Channel.fromPath(params.idx +'/geneInfo.tab')
	  ref_8  = Channel.fromPath(params.idx +'/Genome')
	  ref_9  = Channel.fromPath(params.idx +'/genomeParameters.txt')
	  ref_10 = Channel.fromPath(params.idx +'/SA')
	  ref_11 = Channel.fromPath(params.idx +'/SAindex')
	  ref_12 = Channel.fromPath(params.idx +'/sjdbInfo.txt')
	  ref_13 = Channel.fromPath(params.idx +'/transcriptInfo.tab')
	  ref_14 = Channel.fromPath(params.idx +'/sjdbList.fromGTF.out.tab')
	  ref_15 = Channel.fromPath(params.idx +'/sjdbList.out.tab')
    
    genome_indexed_input = ref_1.concat(ref_2,ref_3,ref_4,ref_5,ref_6,ref_7,ref_8,ref_9,ref_10,ref_11,ref_12,ref_13,ref_14,ref_15)
    // Channel
    //   .of( idx_genome )
    //   .set { genome_indexed_input }
    mapping_withindex(genome_indexed_input.collect(), fastp.out.fastq)
    stats_bam(mapping_withindex.out.bam)
    filter_bam_mapped(mapping_withindex.out.bam)
  }

  sort_bam(filter_bam_mapped.out.bam)
  index_bam(sort_bam.out.bam)
  htseq_count(index_bam.out.bam_idx, gtf_file)

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
}
