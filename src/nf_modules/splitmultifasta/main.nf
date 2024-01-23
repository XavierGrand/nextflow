version = "1.0"
container_url = "xgrand/splitmultifasta:${version}"

params.splitmultifasta_out = ""
process splitmultifasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${fasta.baseName}"
  if (params.splitmultifasta_out != "") {
    publishDir "results/${params.splitmultifasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    path("splited/*.fasta"), emit: splitedfasta
    
  script:
    """
    mkdir splited
    python /app/splitmultifasta.py ${fasta} splited/
    """
}

params.groupsfasta_out = ""
process groupsfasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${genotype}_Genomes"
  if (params.groupsfasta_out != "") {
    publishDir "results/${params.groupsfasta_out}", mode: 'copy'
  }

  input:
    tuple val(genotype), path(fastafiles)

  output:
    tuple val(genotype), path("${genotype}_Genomes_doubled.fasta"), emit: groupedfasta
    
  script:
    """
    touch ${genotype}_Genomes_doubled.fasta
    cat ${fastafiles} >> ${genotype}_Genomes_doubled.fasta
    """
}