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


params.fastp = ""
if(params.fastp == "") {
    if (!params.paired_end) {
        params_fastp = "-x 10 -3 --cut_tail_window_size 10 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence=AAAAAA"
    } else {
        params_fastp = "-x 10 -3 --cut_tail_window_size 10 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence=AAAAAA -F 18 --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    }
} else {
    params_fastp = params.fastp
}
/* Fastp parameters used to customize trimming and dapter removal

@ Type: String
*/


params.spikein_fasta = ""
/* An optional argument specifing a fasta file containing spike-in sequences

Both of those argument must be left empty or filled to use
@Type: String
*/
params.spikein_gtf = ""
/* An optional argument specifing a gtf file for spike-in sequences

@Type: String
*/


if (params.spikein_fasta != "" &&  params.spikein_gtf == "" || 
    params.spikein_fasta == "" &&  params.spikein_gtf != "" ) {
    println "\033[91mError: Only one of the parameters spikein_fasta and \
spikein_gtf were set !\u001B[0m"
    System.exit(1)
} 
spike_in_analysis = params.spikein_fasta != "" && params.spikein_gtf != ""

/*
 ****************************************************************
                              Logs
 ****************************************************************
*/


log.info "paired-end data: ${params.paired_end}"
log.info "fastq files : ${params.fastq}"
log.info "genome files : ${params.genome}"
log.info "htseq param: -s ${params.htseq_param}"
log.info "gtf file: ${params.gtf}"
log.info "index file: ${params.index}"
log.info "fastp parameters: ${params_fastp}"
log.info "spike-in fasta: ${params.spikein_fasta}"
log.info "spike-in gtf: ${params.spikein_gtf}"


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

if (spike_in_analysis) {
    Channel
        .fromPath( params.spikein_fasta )
        .ifEmpty { error "Cannot load spike-in fasta files matching : ${params.spikein_fasta}"}
        .map( it -> [it.baseName, it])
        .set { spikein_fasta_file }
    
    Channel
        .fromPath( params.spikein_gtf )
        .ifEmpty { error "Cannot load spike-in gtf files matching : ${params.spikein_gtf}"}
        .set {spikein_gtf_file }


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


/*   ****************************** Main imports  ******************************************* */

include { fastp_default } from fastp_mod addParams(fastp_out: '02_fastp', 
                                                   fastp: "${params_fastp}" )
include { genome_mapping; index_fasta } from hisat2_mod addParams(folder: '06_mapping')
include { sort_bam } from sammod addParams(sort_bam_out: '07_sort_bam' )
include { index_bam } from sammod addParams(index_bam_out: '08_index_bam' )
include { htseq_count } from htseq_mod addParams(htseq_out: '09_htseq_count', htseq_param: "${params.htseq_param}" )


include { fastqc_fastq as fastqc1} from fastqc_mod addParams(fastqc_fastq_out: '01_fastqc')
include { fastqc_fastq as fastqc2} from fastqc_mod addParams(fastqc_fastq_out: '03_fastqc_trimmed')
include { fastqc_fastq as fastqc3} from fastqc_mod addParams(fastqc_fastq_out: '04_fastqc_nospikein')
include { fastqc_fastq as fastqc_aligned } from fastqc_mod addParams(fastqc_fastq_out: '06_fastqc_aligned')
include { fastqc_fastq as fastqc_unaligned} from fastqc_mod addParams(fastqc_fastq_out: '06_fastqc_unaligned')
include { multiqc_default } from multiqc_mod addParams(multiqc: '--interactive', multiqc_out: "06_multiqc")

/*   ****************************** Spike-in imports  ******************************************* */

include { genome_mapping as genome_mapping_spike; index_fasta as index_fasta_spike } from hisat2_mod addParams(folder: 'spike-in/01_mapping_spike-in',
                                                                                                               notaligned_name: 'no_spike-in')
if (params.paired_end) {
    include { bam_to_fastq_pairedend as bam_to_fastq_spike} from  bedt_mod addParams(bam_to_fastq_pairedend_out: 'spike-in/02_bam_to_fastq_out')
    include { sort_bam as sort_bam_pe_spike } from sammod addParams(sort_bam_out: 'spike-in/02_bamsortname', sort_bam: "-n" )
} else {
    include { bam_to_fastq_singleend as bam_to_fastq_spike} from  bedt_mod addParams(bam_to_fastq_singleend_out: 'spike-in/02_bam_to_fastq_out')
}
include { sort_bam as sort_bam_spike } from sammod addParams(sort_bam_out: 'spike-in/03_sort_bam' )
include { index_bam as index_bam_spike } from sammod addParams(index_bam_out: 'spike-in/04_index_bam' )
include { htseq_count as htseq_count_spike } from htseq_mod addParams(htseq_out: 'spike-in/05_htseq_count', htseq_param: "${params.htseq_param}" )

// spike-in qc
include { fastqc_fastq as fastqc_aligned_spike } from fastqc_mod addParams(fastqc_fastq_out: 'spike-in/01_fastqc_aligned')
include { fastqc_fastq as fastqc_unaligned_spike } from fastqc_mod addParams(fastqc_fastq_out: 'spike-in/01_fastqc_unaligned')
include { multiqc_default as multiqc_default_spike } from multiqc_mod addParams(multiqc: '--interactive',
                                                             multiqc_out: "spike-in/01_multiqc")


/*
 ****************************************************************
                          Workflow
 ****************************************************************
*/

def merge_channels(report1, report2, report3, report4, 
                   fastp_report, hisat2_report) {
    return report1.concat(report2, report3, report4,
    fastp_report.map { it -> [it[0], [it[1]]]},
    hisat2_report.map { it -> [it.baseName, [it]]}).map {it -> it[1]}.collect()
}

workflow spikein_analysis {
    take:
        gtf_spike_in
        genome_spike_in
        fastq_spike_in
    main:
        index_fasta_spike(genome_spike_in)
        genome_mapping_spike(fastq_spike_in, index_fasta_spike.out.index.collect())
        if (!params.paired_end) {
            bam_to_fastq_spike(genome_mapping_spike.out.aligned)
        } else {
            sort_bam_pe_spike(genome_mapping_spike.out.aligned)
            bam_to_fastq_spike(sort_bam_spike.out.bam)
        }
        sort_bam_spike(genome_mapping_spike.out.aligned)
        index_bam_spike(sort_bam_spike.out.bam)
        htseq_count_spike(index_bam_spike.out.bam_idx, gtf_spike_in.collect())

        // QC
        fastqc_unaligned_spike(genome_mapping_spike.out.unaligned)
        fastqc_aligned_spike(bam_to_fastq_spike.out.fastq)
        res = fastqc_unaligned_spike.out.report
            .concat(fastqc_aligned_spike.out.report)
            .map {it -> it[1]}
            .collect()
        multiqc_default_spike(res)
    emit:
        fastq_filtered = genome_mapping_spike.out.unaligned
}


workflow {
    fastp_default(fastq_files)
    if ( params.index == "" ) {
        index_fasta(genome_file)
        index_fasta.out.index.set { index_file }
    }
    if (params.spikein_fasta != "" && params.spikein_gtf != "") {
        spikein_analysis(spikein_gtf_file, spikein_fasta_file, fastp_default.out.fastq)
        spikein_analysis.out.fastq_filtered.map{ it -> [it[1].simpleName, it[1]]}.set { fastq_mapping }
        fastqc3(fastq_mapping)
        fastqc3.out.report.set { fastqc_spikein_report }
    } else {
        fastp_default.out.fastq.set { fastq_mapping }
        Channel
            .fromPath( "tmp_void.txt" )
            .map( it -> [ "tmp", [it]] )
            .set { fastqc_spikein_report }
    }
    genome_mapping(fastq_mapping, index_file.collect())
    sort_bam(genome_mapping.out.aligned)
    index_bam(sort_bam.out.bam)
    htseq_count(index_bam.out.bam_idx, gtf_file.collect())

    // Quality control
    fastqc1(fastq_files)
    fastqc2(fastp_default.out.fastq)
    fastqc_unaligned(genome_mapping.out.unaligned)

    // multiqc
    res = merge_channels(fastqc1.out.report, fastqc2.out.report, fastqc_unaligned.out.report,
                         fastqc_spikein_report, fastp_default.out.report,
                         genome_mapping.out.report)
    multiqc_default(res)
}
