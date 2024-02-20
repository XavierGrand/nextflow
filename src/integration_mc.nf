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
 
/* Define some arguments for cli */
params.project = ""
params.fastq = "${params.project}/fastq/*_{R1,R2}_*.fq.gz"
params.fasta = ""
params.idx = ""

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

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
include { index_fasta as index_fasta_bwa  } from "./nf_modules/bwa/0.7.17/main.nf"
include { mapping_fastq as mapping_fastq_bwa  } from "./nf_modules/bwa/0.7.17/main.nf" addParams(mapping_fastq_out: "01_original_alignment")
include { mapping_fastq as mapping_fastq_softclip  } from "./nf_modules/bwa/0.7.17/main.nf" addParams(mapping_fastq_out: "02_soft_clip_alignment", file_suffix: "_softclip")
include { index_bam } from "./nf_modules/sambamba/1.0.1/main.nf"
include { mark_dup } from "./nf_modules/sambamba/1.0.1/main.nf"
include { sort_bam } from "./nf_modules/sambamba/1.0.1/main.nf"
include { get_soft_clipped } from "./nf_modules/seekSV/1.2.3/main.nf"
include { get_sv } from "./nf_modules/seekSV/1.2.3/main.nf" addParams(get_sv_out: "03_sv_detection")
/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  //############ GENOME INDEXING AND MAPPING ###################

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
      .ifEmpty { error "Cannot find any files matching: ${params.idx}" }
      .map{it -> [(it.baseName =~ /([^\.]*)/)[0][1],[it]]}
      .set { genome_indexed_input }

    mapping_fastq_bwa(genome_indexed_input, fastq_files)
  }

  //#####################DUPLICATE MARKING
  mark_dup(mapping_fastq_bwa.out.bam)
  //#####################COORDINATE SORTING
  sort_bam(mark_dup.out.bam)
  //#####################BAM INDEXING
  //this is not needed since the sorting from sambamba generates already the index
  //index_bam(sort_bam.out.bam)
  //#####################SOFT CLIPPED READS EXTRACTION
/*  Channel
    .fromFilePairs(sort_bam.out.bam)
    .set(sorted_bam_out)*/

  //sort_bam_out=sort_bam.out.bam.map{it.first()}
  //sort_bam.out.view()
  //get_soft_clipped(sort_bam.out.bam)
  get_soft_clipped(sort_bam.out.bam)
  //get_soft_clipped.out.clip_fq.view()
  //index_fasta_bwa.out.index.view()
  //get_soft_clipped(sort_bam_out)
  //#####################REALIGNMENT
  //Need to use the output of index_fasta_bwa or pass the fasta and indexed files as parameters
  // get only the faasta with the clipped reads
  if (params.idx == "") {
   /* Channel
      .fromPath( params.fasta )
      .ifEmpty { error "Cannot find any files matching: ${params.fasta}" }
      .map{it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
      .set { genome_file }
   */ 
    mapping_fastq_softclip(index_fasta_bwa.out.index.collect(), get_soft_clipped.out.clip_fq)
  }
  else {
    //this need to be fixed to mimic the output of the index process
    idx_genome = "${params.idx}"
    
    Channel
      .fromPath( "${params.idx}" )
      .set { genome_indexed_input }

    mapping_fastq_softclip(genome_indexed_input.collect(), get_soft_clipped.out.clip_fq)
  }

  //#####################GET SV
  get_sv(mapping_fastq_softclip.out.bam,sort_bam.out.bam,get_soft_clipped.out.clip_gz)
}
