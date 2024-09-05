#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================================================
                                                      quantseq2
========================================================================================================================

star_fusion pipeline :
 * Pipeline dedicated to transcriptomic data analysis from short-reads.

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
      --design [path]                 Path to the design file.

    References:                       Can be downloaded with download_references.sh (not implemented in pipeline).
      --genome [path]                 Path to genome reference fasta file.
      --index [path]                  Path to STAR-indexed genome folder.
      --gtf [path]                    Path to genome annotation gtf file.

    Optionnal:
      --cyto [path]                   Path to the cytobands file. To draw circos plot by Arriba draw_fusion.R.
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
 
params.fastq = "./data/fastq/*_R{1,2}.fastq.gz"
params.bam = ""
params.design = ""
params.genome = ""
params.cyto = ""
params.gtf = ""
params.index = ""
params.htseq_param = "yes"

/*
 ****************************************************************
                      DESEQ2 PARAMETERS
 ****************************************************************
*/

params.design = ""
/* The file containing the design to use for deseq2 differential expression
analysis.
See https://gitbio.ens-lyon.fr/LBMC/regards/deseq2-wrapper for details
This parameter is optional. If not set: the differential expression process is skipped

@type: File
*/

params.filter = ""
/* An optional parameter used to only keep a subset of gene for differential
expression analysis
See https://gitbio.ens-lyon.fr/LBMC/regards/deseq2-wrapper for details
@type: Optional File
*/

params.gene_name = ""
/* An optional file that can be used to link gene id (ex ensembl! gene id) to
gene symbol
See https://gitbio.ens-lyon.fr/LBMC/regards/deseq2-wrapper for details
@type: Optional file
*/

params.formula = "~ condition"
/*
The formula to use for differential expression, the columns specified here
should be present in the design file. See
https://gitbio.ens-lyon.fr/LBMC/regards/deseq2-wrapper for details

@type: String
*/

params.comparison = "condition-TEST-CTRL"
/* The comparisons to perform. See
https://gitbio.ens-lyon.fr/LBMC/regards/deseq2-wrapper for details

@type: String
*/

params.gene_expression_threshold = 2
/* Minimum expression a gene should have to be concidered prior the diffrential
expression analysis

@type: int
*/

params.basemean_threshold = 0
/*   A threshold to keep significant genes with at least this baseMean
(after deseq2 analysis)
@type: int
*/

params.lfc_threshold = 0
/*
 A threshold to keep significant genes with at least this log2fc
@type: int
*/

/* Specific Arriba parameters
params.arriba_options = "-f 'blacklist,read_through' \
       -E 1 \
       -R 5000 \
       -A 15 \
       -M 2 \
       -U 1000"
*/

/* Params out */
params.fastp_out = "02_fastp"
params.star_index_out = "05_Indexed_genome"
params.star_mapping_fastq_out = "06_mapping_fastq"
params.filter_bam_quality_out = "07_Filtered_bam"
params.index_bam_out = "07_Filtered_bam"
params.deseq2_out = "08_deseq2"
params.star_mapping2fusion_out = "09_mapping2fusion"
params.arriba_out = "10_Arriba_results"
params.draw_fusions_out = "11_drawn_fusions"
params.concat_fusion_out = "10_Arriba_results"
params.fusion_out = "13_DFG"

if (params.design != "") {
    log.info "\n**** DESEQ2 parameter ****"

    log.info "Design file: ${params.design}"
    log.info "Optional file used to filter genes in DE analysis : ${params.filter}"
    log.info "Optional file linking gene id to symbol : ${params.gene_name}"
    log.info "Formula for DE ${params.formula}"
    log.info "Comparisons for DE ${params.comparison}"
    log.info "Gene expresion threshold filter (prior DE): ${params.gene_expression_threshold}"
    log.info "Gene expresion threshold filter (after DE): ${params.basemean_threshold}"
    log.info "Log2foldchange threshold: ${params.lfc_threshold}"
}

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
} else {
  fastq_list = "${params.fastq}"
  log.info "Loaded fastq files (--fastq): ${fastq_list}"
}
if(params.index) {
  index_list = "${params.index}/*"
} else {
  index_list = ""
}
log.info "htseq param: -s ${params.htseq_param}"

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
    .set { gtf_file }

if ( params.design ) {
  Channel
    .fromPath( params.design )
    .set { design }
}

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

/* index file */
if (params.index != "") {
    Channel
        .fromPath( index_list )
        .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
        .collect()
        .map { it -> [ "STAR_index", it ]}
        .set { index_file }
}

if (params.design != "") {

    Channel
        .fromPath( params.design )
        .ifEmpty { error "Cannot find design file matching : ${params.design}"}
        .set { design_file }

    pfilter = params.filter == "" ? "NONE" : params.filter
    Channel
        .fromPath( pfilter )
        .set { filter_file }

    pgene_name = params.gene_name == "" ? "EMPTY" : params.gene_name
    Channel
        .fromPath( pgene_name )
        .set { gene_name_file }

}

/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

include { fastp } from './nf_modules/fastp/main.nf'
include { fastqc_fastq_multi as fastqc_raw } from './nf_modules/fastqc/main.nf' addParams(fastqc_fastq_out: '01_fastqc')
include { fastqc_fastq_multi as fastqc_preprocessed } from './nf_modules/fastqc/main.nf' addParams(fastqc_fastq_out: '03_fastqc_trimmed')
include { multiqc } from './nf_modules/multiqc/main.nf' addParams(multiqc_out: '04_multiqc')
include { filter_bam_quality } from './nf_modules/samtools/1.20/main.nf'
include { index_with_gtf } from './nf_modules/star/2.7.11b--h43eeafb_2/main.nf'
include { mapping_fastq } from './nf_modules/star/2.7.11b--h43eeafb_2/main.nf'
include { mapping2fusion } from './nf_modules/star/2.7.11b--h43eeafb_2/main.nf'
include { index_bam } from './nf_modules/samtools/1.20/main.nf'
include { htseq_count } from './nf_modules/htseq/main.nf' addParams(htseq_out: '09_htseq_count', htseq_param: "${params.htseq_param}" )
include { arriba } from "./nf_modules/arriba/main.nf"
include { draw_fusions } from "./nf_modules/arriba/main.nf"
include { concat_fusion } from "./nf_modules/concatenate/main.nf"
include { parsefusion } from "./nf_modules/fusion_parser/main.nf"
if (params.design != "") {
    include { deseq2_analysis } from "./nf_modules/deseq2/main.nf" addParams(
        formula: "${params.formula}", comparison: "${params.comparison}",
        gene_expression_threshold: "${params.gene_expression_threshold}",
        basemean_threshold: "${params.basemean_threshold}",
        lfc_threshold: "${params.lfc_threshold}")
}

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

workflow differential_expression {
    take:
        design_file
        filter_file
        gene_name_file
        counts
    main:
    if (params.design != "") {
        deseq2_analysis(design_file, filter_file.collect(),
        gene_name_file.collect(), counts.collect())
        if (params.top > 0) {
            topgo_analysis(deseq2_analysis.out.results)
        }

    }
}

def merge_channels(report1, report2) {
    return report1.concat(report2)
}

workflow {

/*
 ****************************************************************
                      Preprocessing & QC
 ****************************************************************
*/

  if(params.fastq != ""){
    fastp(fastq_files)
    fastqc_raw(fastq_files)
    fastqc_preprocessed(fastp.out.fastq)
    res = merge_channels(fastqc_raw.out.report, fastqc_preprocessed.out.report)
    multiqc(res)

/*
 ****************************************************************
                          Mapping
 ****************************************************************
*/

    if(params.index == "") {
      index_with_gtf(genome_file, gtf_file)
      mapping_fastq(index_with_gtf.out.index.collect(), fastp.out.fastq)
    } else {
      mapping_fastq(index_file.collect(), fastp.out.fastq)
    }

/*
 ****************************************************************
                    Filtering & indexing
 ****************************************************************
*/

    filter_bam_quality(mapping_fastq.out.bam)
    index_bam(filter_bam_quality.out.bam)
  }

  else { index_bam(bam_files) }

/*
 ****************************************************************
                              Count
 ****************************************************************
*/

  htseq_count(index_bam.out.bam_idx, gtf_file.collect())

/*
 ****************************************************************
                Differential Expression analysis
 ****************************************************************
*/

  if (params.design != "") {
    differential_expression(design_file, filter_file, gene_name_file, htseq_count.out.counts)
  }


/*

 ****************************************************************
                      Fusion detection
 ****************************************************************


  arriba(index_bam.out.bam_idx, gtf_file.collect(), genome_file.collect())
  if (params.cyto != "") { draw_fusions(arriba.out.fusions, index_bam.out.bam_idx, gtf_file.collect(), cytobands.collect()) }


 ****************************************************************
              Parsing results & statistical analysis
 ****************************************************************

  if ( params.design ) {
    concat_fusion(arriba.out.fusions)
    parsefusion(concat_fusion.out.concatenated_fusions.collect(), htseq_count.out.counts.collect(), design)
  }

*/

}