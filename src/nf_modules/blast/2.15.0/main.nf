version = "2.15.0"
container_url = "ncbi/blast:${version}"

params.makeblastdb_out = ""
process makeblastdb {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}_Genomes"
  if (params.makeblastdb_out != "") {
    publishDir "results/${params.makeblastdb_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(ref_fasta)

  output:
    tuple val(ref_fasta.baseName), path("*.fasta.n*"), emit: blastdb

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

params.blast_them_all_out = ""
process blast_them_all {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${file_id}_Genomes"
  if (params.blast_them_all_out != "") {
    publishDir "results/${params.blast_them_all_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(genotype), path(blastdb)

  output:
    path("${file_id}_hits.txt"), emit: blasthits
    path("${file_id}_hits_counts.txt"), emit: counthits
    path("${file_id}_best_ref.txt"), emit: bestref

  script:
"""
blastn -db ${genotype}.fasta -query ${fastq} \
       -task megablast \
       -max_target_seqs 1 \
       -max_hsps 1 \
			 -outfmt "6 qseqid sseqid evalue bitscore slen qlen length pident" \
			 -out ${file_id}_hits.txt -num_threads ${params.blasthreads}
cut -f2 ${file_id}_hits.txt | sort | uniq -c | sort -k 1,1 -r > ${file_id}_hits_counts.txt
cut -f2 ${file_id}_hits.txt | sort | uniq -c | sort -k 1,1 -r | head -n1 | sed 's/^ *[0-9]* //g' > ${file_id}_best_ref.txt
"""
}
