nextflow.enable.dsl=2

/*
Testing pipeline for marseq scRNASeq analysis
*/

include { adaptor_removal} from "./nf_modules/cutadapt/main.nf"
include {
  index_fasta;
  count;
  index_fasta_velocity;
  count_velocity
} from "./nf_modules/kb/main.nf" addParams(
  kb_protocol: "marsseq",
  count_out: "quantification/",
  count_velocity_out: "quantification_velocity/"
)

params.fasta = "http://ftp.ensembl.org/pub/release-98/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz"
params.fastq = "data/CF42_45/*/*R{1,2}.fastq.gz"
params.gtf = "http://ftp.ensembl.org/pub/release-98/gtf/gallus_gallus/Gallus_gallus.GRCg6a.98.gtf.gz"
params.transcript_to_gene = ""
params.whitelist = "data/expected_whitelist.txt"
params.config = "data/marseq_flexi_splitter.yaml"
params.workflow_type = "classic"

log.info "fastq files (--fastq): ${params.fastq}"
log.info "fasta file (--fasta): ${params.fasta}"
log.info "gtf file (--gtf): ${params.gtf}"
log.info "transcript_to_gene file (--transcript_to_gene): ${params.transcript_to_gene}"
log.info "whitelist file (--whitelist): ${params.whitelist}"
log.info "config file (--config): ${params.config}"

channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }
channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }
channel
  .fromPath( params.gtf )
  .ifEmpty { error "Cannot find any gtf files matching: ${params.gtf}" }
  .map { it -> [it.simpleName, it]}
  .set { gtf_files }
if (params.whitelist == "") {
  channel.empty()
    .set { whitelist_files }
} else {
  channel
    .fromPath( params.whitelist )
    .map { it -> [it.simpleName, it]}
    .set { whitelist_files }
}
channel
  .fromPath( params.config )
  .ifEmpty { error "Cannot find any config files matching: ${params.config}" }
  .map { it -> [it.simpleName, it]}
  .set { config_files }

workflow {
  adaptor_removal(fastq_files)
  if (params.workflow_type == "classic") {
    index_fasta(
      fasta_files,
      gtf_files
    )
    count(
      index_fasta.out.index,
      adaptor_removal.out.fastq,
      index_fasta.out.t2g, whitelist_files,config_files
    )
  } else {
    index_fasta_velocity(
      fasta_files,
      gtf_files
    )
    count_velocity(
      index_fasta_velocity.out.index,
      adaptor_removal.out.fastq,
      index_fasta_velocity.out.t2g,
      whitelist_files,
      config_files
    )
  }
}
