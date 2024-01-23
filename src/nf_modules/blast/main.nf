version = "2.15.0"
container_url = "ncbi/blast:${version}"

params.ref_fasta = ""
params.makeblastdb_out = ""
process makeblastdb {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.makeblastdb_out != "") {
    publishDir "results/${params.makeblastdb_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(ref_fasta)

  output:
    tuple val(file_id), path("*.fasta.n*"), emit: blastdb

  script:
"""
makeblastdb \
    -in ${ref_fasta} \
    -input_type 'fasta' \
    -dbtype 'nucl' \
    -parse_seqids
"""
}

params.hbvdb = ""
process dl_hbvdb {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "${params.hbvdb}_Genomes"
  if (params.dl_hbvdb_out != "") {
    publishDir "results/${params.dl_hbvdb_out}", mode: 'copy'
  }

  input:

  output:
    tuple val(params.hbvdb), path("*.fasta"), emit: reference_db

  script:
  if (params.hbvdb == "all" || params.hbvdb == "A" || params.hbvdb == "B" || params.hbvdb == "C" || 
      params.hbvdb == "D" || params.hbvdb == "E" || params.hbvdb == "F" || params.hbvdb == "G" || 
      params.hbvdb == "H") { // || params.hbvdb == "I" || params.hbvdb == "J") {
    link = "https://hbvdb.lyon.inserm.fr/data/nucleic/fasta/${params.hbvdb}_Genomes.fas"
    output_name = "${params.hbvdb}_Genomes.fasta"
  }
  // else if(){}
"""
wget --quiet --no-check-certificate -O ${output_name} ${link}
"""
}