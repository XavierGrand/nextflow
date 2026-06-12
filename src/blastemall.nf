#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
========================================================
                    Blast'em All
========================================================

Blast'em All pipeline :
  - Blast any sequence against any database

********************************************************
                    Help message definition
********************************************************
*/

def help() {
  log.info"""
  Usage:
  The minimal command for running the pipeline is as follows:

    nextflow src/majora.nf -c src/nextflow.config -profile <profile> --fastq or --fast5 or --pod5 <path/to/files> --ref_db <path/to/csv> or ---hbvdb [string] or --ref_user <path/to/fasta_or_multifasta>

  Nextflow parameters:
    -profile [str]                  Configuration profile to use. MANDATORY
                                    Available: docker, singularity, podman, psmn, ccin2p3, etc.
                                    User can use his own, associated with a nextflow configuration file.
                                    Refer to https://www.nextflow.io/docs/latest/config.html#configuration-file. 

  Input file:                       Only one of the three input types can be specified. MANDATORY
    --query [path]                  Path to the query files.
    
  References:                       These sequences will be used during BLAST alignment. MANDATORY
  --ref_db [path]                   To provide csv (tab separated) file containg hbvdb sequence IDs in the second column.
                                      
  --hbvdb  [str]                    Should be followed by either "all" or the uppercase letter corresponding to the desired genotype.
                                      
  --ref_user [path]                 To provide fasta/multifasta file containg hbvdb sequences.

  Help:
    --help | --h                    Display this help message.
    
  """.stripIndent()
}

// Show help message
params.help = '' // initialisation du paramètre help à une chaine 
params.h = ''

if (params.help || params.h) {
  help()
  exit 0
}

/*
********************************************************
                   Default Parameters
********************************************************
*/

/*
*** Params in ***
*/

params.blasthreads = 18
params.query = ""
params.hbvdb = ""
params.ref_db = ""
params.ref_user = ""

params.dl_hbvdb_out =  '01_dl_hdvdb/'
params.makerefdb_out = '02_makerefdb/'
params.makeblastdb_out = '03_makeblastdb/'
params.blast_them_all_out = '04_blast_them_all/'

/*
 ****************************************************************
                        Channel Definitions
 ****************************************************************
*/

// Create input channel depending on input file extension

// QUERY input
if (params.query != '') {
  Channel.fromPath(params.query)
         .map(it -> [it.baseName, it])
         .set { query }
}

// To download hbv or hdv data bases (all or selected genome)
if (params.hbvdb != '') {
  Channel.fromPath(params.hbvdb)
         .set { hbvdb }
}


// To make the user use its own multi-fasta file
if (params.ref_user != '') {
  Channel.fromPath(params.ref_user)
         .map(it -> [it.baseName, it])
         .set { user_ref }
}


// To use a csv file containing reference sequences ID from hbvdb
if (params.ref_db != '') {
  Channel.fromPath(params.ref_db)
         .set { ref_db }
}

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { dl_hbvdb } from './nf_modules/blast/2.15.0/main.nf'
include { makeblastdb } from './nf_modules/blast/2.15.0/main.nf'
include { makerefdb } from './nf_modules/blast/2.15.0/main.nf'
include { blast_them_all } from './nf_modules/blast/2.15.0/main.nf'

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/
workflow {
/*
**********************Download/upload reference depending on used option *********************
*/

  if ( params.ref_db != '' ) {
      makerefdb(ref_db)
      makeblastdb(makerefdb.out.ref_db)
  }

  else if ( params.hbvdb != '' ) {
    dl_hbvdb()
    dl_hbvdb.out.reference_db.set { genotype }
    makeblastdb(genotype)
  }

  else if ( params.ref_user != '' ) {
   makeblastdb(user_ref)
  }

  blast_them_all(query, makeblastdb.out.blastdb.collect())

}