#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
##################################################################
Integration pipeline
##################################################################

Pipeline to perform viral integration analysis and discovery, using both DNA and RNA data.


/*
****************************************************************
                    Default Parameters
****************************************************************
*/
 
/* Arguments */
params.project = ""
params.fastq = "${params.project}/fastq/*_{R1,R2}_*.fq.gz"
params.gtf = ""
params.fasta = ""
params.idx = ""
params.filter_bam_mapped = "-F 268 -f 1 -q 10"

params.star_mapping_fastq_out = "STAR/"
params.star_index_out = "STARindex/"
params.sort_bam_out = "STAR/"
params.filter_bam_out = "STAR/"

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


/*
 ****************************************************************
                          Imports
 ****************************************************************
*/
include { mapping } from "./nf_modules/bwa/main.nf"
include { index_fasta as index_fasta_bwa  } from "./nf_modules/bwa/main.nf"
include { mapping_fastq as mapping_fastq_bwa  } from "./nf_modules/bwa/main.nf"
/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  //############ GENOME INDEXATION AND MAPPING ###################

  if (params.idx == "") {
    Channel
      .fromPath( params.fasta )
      .ifEmpty { error "Cannot find any files matching: ${params.fasta}" }
      .map{it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
      .set { genome_file }
    
    index_fasta_bwa(genome_file)
    mapping_fastq_bwa(index_fasta_bwa.out.index.collect(), fastq_files)
  }
  else {
    idx_genome = "${params.idx}"
    
    Channel
      .fromPath( "${params.idx}" )
      .set { genome_indexed_input }

    mapping_fastq_bwa(genome_indexed_input, fastq_files)
  }


}
