nextflow.enable.dsl=2

/*
Testing pipeline for marseq scRNASeq analysis
*/

include {
  mapping;
} from "./nf_modules/bwa/main.nf"

include {
  sort_bam;
} from "./nf_modules/samtools/main.nf"

include {
  germline_cohort_data_variant_calling;
} from "./nf_modules/gatk4/main.nf" addParams(
  variant_calling_out: "vcf/",
)

params.fastq = ""
params.fasta = ""

channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }
channel
  .fromPath( params.fasta )
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }

workflow {
  mapping(fasta_files, fastq_files)
  sort_bam(mapping.out.bam)
  germline_cohort_data_variant_calling(sort_bam.out.bam, fasta_files)
}
