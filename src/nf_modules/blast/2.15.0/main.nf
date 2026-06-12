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
    //tuple val(params.hbvdb), path("*.fasta"), emit: reference_db
    path("*.fasta"), emit: reference_db

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
  tag "${barcode}_Genomes"
  if (params.blast_them_all_out != "") {
    publishDir "results/${params.blast_them_all_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(fastq)
    tuple val(genotype), path(blastdb)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_hits.txt"), emit: blasthits
    tuple val(barcode), path("${barcode}/${barcode}_hits_counts.txt"), emit: counthits
    tuple val(barcode), path("${barcode}/${barcode}_best_ref.txt"), emit: bestref

  script:
"""
mkdir ${barcode}
blastn -db ${genotype}.fasta -query ${fastq} \
       -task megablast \
       -max_target_seqs 1 \
       -max_hsps 1 \
			 -outfmt "6 qseqid sseqid evalue bitscore slen qlen length pident" \
			 -out ${barcode}/${barcode}_hits.txt -num_threads ${task.cpus}
cut -f2 ${barcode}/${barcode}_hits.txt | sort | uniq -c | sort -k 1,1 -r > ${barcode}/${barcode}_hits_counts.txt
cut -f2 ${barcode}/${barcode}_hits.txt | sort | uniq -c | sort -k 1,1 -r | head -n1 | sed 's/^ *[0-9]* //g' > ${barcode}/${barcode}_best_ref.txt
"""
}

process makerefdb {
  label "big_mem_multi_cpus"
  // mettre en tag l'ID de la séquence_Genomes
  if (params.makerefdb_out != "") {
    publishDir "results/${params.makerefdb_out}", mode: 'copy'
  }

  input:
  path(ref_file)

  output: 
  path("combined_ref.fasta"), emit: ref_db
  
  script:
"""
  #!/usr/bin/bash 
  
  awk 'BEGIN { FS = "\t" } NR == 1 { for (i = 1; i <= NF; i++) {if (\$i == "ID") { id_col = i}} if (!id_col) { print "The column ID has not been found"} } NR > 1 { print "https://hbvdb.lyon.inserm.fr/tmp/hbvdb_dat/" \$id_col "/" \$id_col "_sequence.txt"}' ${ref_file}  > temporary.txt
  mkdir refdb
  wget -P refdb -i temporary.txt --no-check-certificate
  touch combined_ref.fasta

  for fichier in refdb/*_sequence.txt
  do
    cat "\$fichier" >> combined_ref.fasta
    echo "" >> combined_ref.fasta
  done
  sed -i 's/>/>gnl|hbvnuc|/' combined_ref.fasta
  rm temporary.txt
"""
}
