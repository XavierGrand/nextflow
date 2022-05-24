nextflow.enable.dsl=2

/* Quant-seq pipeline */


/*
 ****************************************************************
                       parameters
 ****************************************************************
*/

params.paired_end = true
/* false for single end data, true for paired-end data

@type: Boolean
*/

params.fastq = "data/tiny_ribosome_profiling_dataset/fastq/*.fastq"
/* Fastq files quant-seq reads

@type: Files
*/

params.genome = "data/genome.fa"
/* A genome file
@type: File
*/

params.htseq_param = "yes"
/* yes/no/reverse whether the data is from a strand-specific assay (default: yes)
@Type: String
*/

params.gtf = "data/test.gtf"
/* A gtf file to use for htseq count
@Type: File
*/


params.index = ""
/* Path leading to a index file

@ Type: file
*/


if (!params.paired_end) {
    params.fastp = "-x 10 -3 --cut_tail_window_size 10 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence=AAAAAA"
} else {
    params.fastp = "-x 10 -3 --cut_tail_window_size 10 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence=AAAAAA -F 18 --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
}


/*
 ****************************************************************
                              Logs
 ****************************************************************
*/


log.info "paired-end data: ${params.paired_end}"
log.info "fastq files : ${params.fastq}"
log.info "genome files : ${params.genome}"


/*
 ****************************************************************
                  Channel definitions
 ****************************************************************
*/


/* Raw reads fastq */
Channel
  .fromFilePairs( params.fastq, size: -1 )
  .set { fastq_files }


/* genome file */
Channel
  .fromPath( params.genome )
  .ifEmpty { error "Cannot find any genome files matching: ${params.genome}" }
  .map( it -> [it.baseName, it])
  .set { genome_file }


/* gtf file */
Channel
    .fromPath( params.gtf )
    .ifEmpty { error "Cannot find any gtf files matching: ${params.gtf}" }
    .set { gtf_file }


/* index file */
if (params.index != "") {
    Channel
        .fromPath( params.index )
        .ifEmpty { error "Cannot find any index files matching: ${params.index}" }
        .collect()
        .map { it -> [ "ht2_index", it ]}
        .set { index_file }
} else {
    Channel.from( "" ).set{ index_file }
}


/*
 ****************************************************************
                          Imports
 ****************************************************************
*/

fastqc_mod = './nf_modules/fastqc/main.nf'
fastp_mod = './nf_modules/fastp/main.nf'
hisat2_mod = './nf_modules/hisat2/main.nf'
multiqc_mod = './nf_modules/multiqc/main.nf'
bedt_mod = './nf_modules/bedtools/main.nf'
htseq_mod = './nf_modules/htseq/main.nf'
sammod = './nf_modules/samtools/main.nf'


include { fastqc_fastq as fastqc1} from fastqc_mod addParams(fastqc_fastq_out: '01_fastqc')
include { fastp_default } from fastp_mod addParams(fastp_out: '02_fastp', 
                                                   fastp: "${params.fastp}" )
include { fastqc_fastq as fastqc2} from fastqc_mod addParams(fastqc_fastq_out: '03_fastqc_trimmed')
include { fastqc_fastq as fastqc_aligned} from fastqc_mod addParams(fastqc_fastq_out: '05_fastqc_aligned')
include { fastqc_fastq as fastqc_unaligned} from fastqc_mod addParams(fastqc_fastq_out: '06_fastqc_unaligned')
if (params.paired_end) {
    include { bam_to_fastq_pairedend as bam_to_fastq} from  bedt_mod addParams(bam_to_fastq_pairedend_out: '05_bam_to_fastq_out')
} else {
    include { bam_to_fastq_singleend as bam_to_fastq} from  bedt_mod addParams(bam_to_fastq_singleend_out: '05_bam_to_fastq_out')
}
include { sort_bam } from sammod addParams(sort_bam_out: '07_sort_bam' )
include { index_bam } from sammod addParams(index_bam_out: '08_index_bam' )
include { htseq_count } from htseq_mod addParams(htseq_out: '09_htseq_count', htseq_param: "${params.htseq_param}" )
include { genome_mapping; index_fasta } from hisat2_mod addParams(folder: '06_mapping')
include { multiqc_default } from multiqc_mod addParams(multiqc: '--interactive')

/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

def merge_channels(report1, report2, report3, report4) {
    return report1.concat(report2, report3, report4).map {it -> it[1]} 
    .collect()
}

workflow {
    fastp_default(fastq_files)
    if ( params.index == "" ) {
        index_fasta(genome_file)
        index_fasta.out.index.set { index_file }
    }
    genome_mapping(fastp_default.out.fastq, index_file.collect())
    bam_to_fastq(genome_mapping.out.aligned)
    sort_bam(genome_mapping.out.aligned)
    index_bam(sort_bam.out.bam)
    htseq_count(index_bam.out.bam_idx, gtf_file.collect())

    // Quality control
    fastqc1(fastq_files)
    fastqc2(fastp_default.out.fastq)
    fastqc_unaligned(genome_mapping.out.unaligned)
    fastqc_aligned(bam_to_fastq.out.fastq)

    // multiqc
    res = merge_channels(fastqc1.out.report, fastqc2.out.report, fastqc_unaligned.out.report,
                         fastqc_aligned.out.report)
    multiqc_default(res)
}