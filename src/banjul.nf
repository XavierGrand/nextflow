#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                Banjul
========================================================================================================================

Banjul analysis pipeline :
 * Pipeline dedicated to HBV infected cohort genotyping analysis from Nanopore sequencing data.

Author Xavier Grand <xavier.grand@inserm.fr>

Creation date: 16/01/2024 by Xavier Grand

 ****************************************************************
                      Help Message Definition
 ****************************************************************
*/

def helpMessage() {
    log.info"""
    Usage:
    Pipeline dedicated to HBV infected cohort genotyping analysis from Nanopore sequencing data.
    The typical command for running the pipeline is as follows:

      nextflow ./src/banjul.nf -c ./src/nextflow.config -profile singularity

    Configuration argument:
      -profile [str]                Configuration profile to use.
                                    Available: docker, singularity, podman, psmn, ccin2p3

    Samples:
      --input [path]                Path and Pattern to fastq files.

    References:
      --hbvdb [string]              Source of reference genomes.
                                    Available: "all" (default),"A","B","C","D","E","F","G","H","I","J".

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

params.input = "/home/xavier/Data/Gambia_cohort/banjul/01_basecalling/"
params.fasta = ""
params.hbvdb = "all"
params.blasthreads = 18 // To override the nextflow.config threads number configuraton

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/

log.info "Input path: ${params.input}"
log.info "Reference genome fasta file: ${params.fasta}"

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
  .fromPath(params.input+'*/', type: 'dir')
  .map(it -> [it.baseName, it])
  .set{barcodes}

/*
Channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map( it -> [it.baseName, it])
  .set { fasta_file }
*/

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { dl_hbvdb } from "./nf_modules/blast/2.15.0/main.nf" addParams(dl_hbvdb_out: "00_ReferenceFiles/")
include { splitmultifasta } from "./nf_modules/splitmultifasta/1.0/main.nf" addParams(splitmultifasta_out: "00_ReferenceFiles/")
include { doublefastaref } from "./nf_modules/seqkit/2.1.0/main.nf" addParams(doublefastaref_out: "00_ReferenceFiles/")
include { groupsfasta } from "./nf_modules/splitmultifasta/1.0/main.nf" addParams(groupsfasta_out: "00_ReferenceFiles/")
include { makeblastdb } from "./nf_modules/blast/2.15.0/main.nf" addParams(makeblastdb_out: "00_ReferenceFiles/")
include { concatenate } from "./nf_modules/seqkit/2.1.0/main.nf" addParams(fastq_out: "01_fastq/")
include { sample_fastq } from "./nf_modules/seqtk/1.3/main.nf" addParams(sample_fastq_out: "02_Blast/")
include { blast_them_all } from "./nf_modules/blast/2.15.0/main.nf" addParams(blast_them_all_out: "02_Blast/")
include { extractref } from "./nf_modules/seqkit/2.1.0/main.nf" addParams(extractref_out: "00_ReferenceFiles/")
include { index_fasta } from "./nf_modules/samtools/1.11/main.nf" addParams(index_fasta_out: "00_ReferenceFiles/")
include { mapping_hbv_genome } from "./nf_modules/minimap2/2.17/main.nf" addParams(mapping_hbv_genome_out: "03_Minimap2/")

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  concatenate(barcodes)
  sample_fastq(concatenate.out.merged_fastq)

// Step 1: Download or Upload Reference to blast:
/*
Download HBV reference sequences from hbvdb:
Use blast image, ~ Ubuntu 22.04.3 LTS based.

Optionnal: one of these options:
#1 DL all reference from hbvdb:
#2 DL selected genotype(s) only: need a parameter --gen [string] e.g. A,B,E as a Channel genotypes
*/
  if (params.hbvdb != "") {
    dl_hbvdb()
  }

/*
#3 Load user's multi-fasta file: need a parameter --references
Channel in or params...

#4 Load user's blastdb files: need a parameter --blastdb
Channel in or params...


// Step 2 : extract single sequence fasta from multifasta files and make blastdb.
if DL ref from hbvdb OR Load user's multi-fasta file : doubled reference sequences individually using "seqkit concat" then, "makeblastdb"
else Load user's blastdb.
*/

  splitmultifasta(dl_hbvdb.out.reference_db)
  doublefastaref(splitmultifasta.out.splitedfasta)
  groupsfasta(doublefastaref.out.doubledfasta.groupTuple())
  makeblastdb(groupsfasta.out.groupedfasta)

// Step 3 : Blast reads from fastq files on blastdb

  blast_them_all(sample_fastq.out.sampled_fastq, makeblastdb.out.blastdb)

// Step 4 : Extract best results and corresponding reference sequence

  extractref(blast_them_all.out.bestref, groupsfasta.out.groupedfasta)
  index_fasta(extractref.out.referenceseq)

// Step 5 : Align reads on reference sequence

  mapping_hbv_genome(concatenate.out.merged_fastq, extractref.out.referenceseq)

/*
// Step 6 : filter mapping results
samtools, ok

// Step 7 : Variation calling

*/
}