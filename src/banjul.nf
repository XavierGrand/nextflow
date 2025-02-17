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
      --fasta [path]                Optionnal: Path to multi-fasta file containing reference sequences.

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

params.input = "/home/xavier/Data/Gambia_cohort/20230214_minION_BS_HighVL_couple_SHAD01/01_basecalling/"
params.fasta = ""
params.hbvdb = "all"
params.blasthreads = 18 // To override the nextflow.config threads number configuraton
params.fwprimer = "CTACTGTTCAAGCCTCCAAGC" // Denomination P3(1857-1877)
params.rwprimer = "CGCAGACCAATTTATGCCTAC" // Denomination P4(1783-1803)
params.maxlength = 3300
params.minlength = 3000
params.filter_bam_mapped = "-f0 -f16"

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

if ( params.fasta != "" ) {
  Channel
    .fromPath( params.fasta )
    .set { fasta_file }
}

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { dl_hbvdb } from "./nf_modules/blast/2.15.0/main.nf" addParams(dl_hbvdb_out: "00_ReferenceFiles/")
include { splitmultifasta } from "./nf_modules/splitmultifasta/1.0/main.nf" addParams(splitmultifasta_out: "00_ReferenceFiles/")
include { doublefastaref } from "./nf_modules/seqkit/2.8.2/main.nf" addParams(doublefastaref_out: "00_ReferenceFiles/")
include { groupsfasta } from "./nf_modules/splitmultifasta/1.0/main.nf" addParams(groupsfasta_out: "00_ReferenceFiles/")
include { makeblastdb } from "./nf_modules/blast/2.15.0/main.nf" addParams(makeblastdb_out: "00_ReferenceFiles/")
include { concatenate } from "./nf_modules/seqkit/2.8.2/main.nf" addParams(fastq_out: "01_fastq/")
include { grep_primer as pick_fw_primer } from "./nf_modules/seqkit/2.8.2/main.nf" addParams(grep_primer_out: "01_fastq/", primerloc: "fw")
include { grep_primer as pick_rw_primer } from "./nf_modules/seqkit/2.8.2/main.nf" addParams(grep_primer_out: "01_fastq/", primerloc: "rw")
include { porechop } from "./nf_modules/porechop/0.2.4/main.nf" addParams(porechop_out: "01_fastq/")
include { filterbylength } from "./nf_modules/seqkit/2.8.2/main.nf" addParams(filterbylength_out: "01_fastq/")
include { sample_fastq } from "./nf_modules/seqtk/1.3/main.nf" addParams(sample_fastq_out: "02_Blast/")
include { blast_them_all } from "./nf_modules/blast/2.15.0/main.nf" addParams(blast_them_all_out: "02_Blast/")
include { extractref } from "./nf_modules/seqkit/2.8.2/main.nf" addParams(extractref_out: "00_ReferenceFiles/")
include { index_fasta } from "./nf_modules/samtools/1.20/main.nf" addParams(index_fasta_out: "00_ReferenceFiles/")
include { mapping_hbv_genome } from "./nf_modules/minimap2/2.17/main.nf" addParams(mapping_hbv_genome_out: "03_Minimap2/")
include { sort_bam } from "./nf_modules/samtools/1.20/main.nf" addParams(sort_bam_out: "03_Minimap2/")
include { index_bam } from "./nf_modules/samtools/1.20/main.nf" addParams(index_bam_out: "03_Minimap2/")
include { filter_bam_mapped } from "./nf_modules/samtools/1.20/main.nf" addParams(filter_bam_mapped_out: "03_Minimap2/")
include { consensus } from "./nf_modules/samtools/1.20/main.nf" addParams(consensus_out: "05_consensus/")

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow {

  concatenate(barcodes)
  pick_fw_primer(concatenate.out.merged_fastq, params.fwprimer)
  pick_rw_primer(pick_fw_primer.out.filtered_fastq, params.rwprimer)
  porechop(pick_rw_primer.out.filtered_fastq)
  filterbylength(porechop.out.porechoped_fastq)
  sample_fastq(filterbylength.out.filtered_fastq)

// Step 1: Download or Upload Reference to blast:
/*
Download HBV reference sequences from hbvdb:
Use blast image, ~ Ubuntu 22.04.3 LTS based.

Optionnal: one of these options:
#1 DL all reference from hbvdb:
#2 DL selected genotype(s) only: need a parameter --gen [string] e.g. A,B,E as a Channel genotypes
*/
  if ( params.fasta == "" ) {
    if ( params.hbvdb != "" ) {
      dl_hbvdb()
      dl_hbvdb.out.reference_db.set { fasta_file }
    }
  }

/*
#3 Load user's multi-fasta file: need a parameter --references
Channel in or params...

// Ok, Channel under condition.

#4 Load user's blastdb files: need a parameter --blastdb
Channel in or params...


// Step 2 : extract single sequence fasta from multifasta files and make blastdb.
if DL ref from hbvdb OR Load user's multi-fasta file : doubled reference sequences individually using "seqkit concat" then, "makeblastdb"
else Load user's blastdb.
*/

  splitmultifasta(fasta_file)
  doublefastaref(splitmultifasta.out.splitedfasta)
  groupsfasta(doublefastaref.out.doubledfasta.groupTuple())
  makeblastdb(groupsfasta.out.groupedfasta)

// Step 3 : Blast reads from fastq files on blastdb

  blast_them_all(sample_fastq.out.sampled_fastq, makeblastdb.out.blastdb.collect())

// Step 4 : Extract best results and corresponding reference sequence

  extractref(blast_them_all.out.bestref, groupsfasta.out.groupedfasta.collect())
  index_fasta(extractref.out.referenceseq)

// Step 5 : Align reads on reference sequence

  mapping_hbv_genome(filterbylength.out.filtered_fastq.combine(extractref.out.referenceseq, by: 0))

// Step 6 : filter mapping results

  filter_bam_mapped(mapping_hbv_genome.out.bam)
  sort_bam(filter_bam_mapped.out.bam)
  index_bam(sort_bam.out.bam)

// Step 7 : Variation calling/Consensus sequence, USE NANOPOLISH

  consensus(sort_bam.out.bam)

// Step 8 : Phylogenetic tree

  //mkinfile()
}
