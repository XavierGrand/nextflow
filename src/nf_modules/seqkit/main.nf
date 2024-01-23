version = "2.1.0"
container_url = "xgrand/seqkit:${version}"

params.reversecomp_out = ""
process reversecomp {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "rev-comp"
  if (params.reversecomp_out != "") {
    publishDir "results/${params.reversecomp_out}", mode: 'copy'
  }
  
  input:
    val(adapt)

  output:
    path("adapt.fasta"), emit: adapt_fst
    path("adaptRC.fasta"), emit: adaptRC_fst

  script:
    """
    echo ">adapt" >> adapt.fasta
    echo ${adapt} >> adapt.fasta
    seqkit seq adapt.fasta -r -p -r -p -t DNA -v > adaptRC.fasta
    """
}

params.seqkit_grep_out = ""
process seqkit_grep {
  container = "${container_url}"
  label "small_mem_multi_cpus"
  tag "${barcode}"
  if (params.seqkit_grep_out != "") {
    publishDir "results/${params.seqkit_grep_out}", mode: 'copy'
  }
  
  input:
    tuple val(barcode), path(fastq)
    val(adapt)
    val(gsp)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_390bp_filtered_5RACE_GSP.fastq.gz"), emit: filtered_fastq
    path("${barcode}/*.csv")
    path("${barcode}/*.txt")
    path("${barcode}/${barcode}_filtered_5RACE.fastq.gz")
    path("${barcode}/${barcode}_filtered_5RACE_GSP.fastq.gz")

  script:
    lgadapt = Math.round(adapt.size().div(10))
    lggsp = Math.round(gsp.size().div(10))
    """
    mkdir ${barcode}
    cd ${barcode}/
    echo "mismatch allowed to 5'RACE adapter:  ${lgadapt}" > mismatch.txt
    echo "mismatch allowed to Gene Specific primer:  ${lggsp}" >> mismatch.txt
    echo ${adapt} > adapt.txt
    echo ${gsp} > gsp.txt
    seqkit grep -i -f adapt.txt -m ${lgadapt} ../${fastq} -o ${barcode}_filtered_5RACE.fastq.gz -j ${task.cpus}
    seqkit grep -i -f gsp.txt -m ${lggsp} ${barcode}_filtered_5RACE.fastq.gz -o ${barcode}_filtered_5RACE_GSP.fastq.gz -j ${task.cpus}
    seqkit seq --min-len 390 --remove-gaps ${barcode}_filtered_5RACE_GSP.fastq.gz -j ${task.cpus} > ${barcode}_390bp_filtered_5RACE_GSP.fastq
    gzip ${barcode}_390bp_filtered_5RACE_GSP.fastq
    seqkit stats ../${fastq} -T -j ${task.cpus} > ${barcode}_seq_stats.csv
    seqkit stats ${barcode}_filtered_5RACE.fastq.gz -T -j ${task.cpus} | tail -n1 >> ${barcode}_seq_stats.csv
    seqkit stats ${barcode}_filtered_5RACE_GSP.fastq.gz -T -j ${task.cpus} | tail -n1 >> ${barcode}_seq_stats.csv
    seqkit stats ${barcode}_390bp_filtered_5RACE_GSP.fastq.gz -T -j ${task.cpus} | tail -n1 >> ${barcode}_seq_stats.csv
    """
}

params.fastq_out = ""
process concatenate {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "${barcode}"
  if (params.fastq_out != "") {
    publishDir "results/${params.fastq_out}", mode: 'copy'
  }

  input:
    tuple val(barcode), path(fastq)

  output:
    tuple val(barcode), path("${barcode}/${barcode}_merged.fastq.gz"), emit: merged_fastq

  script:
    """
    mv ${fastq} path_${fastq}
    mkdir ${barcode}
    cd ${barcode}/
    path=\$(readlink -f ../path_${fastq})
    seqkit scat -j ${task.cpus} -f \${path} --gz-only > ${barcode}_merged.fastq
    gzip ${barcode}_merged.fastq
    """
}

params.doublefastaref_out = ""
process doublefastaref {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "${fasta.baseName}"
  if (params.doublefastaref_out != "") {
    publishDir "results/${params.doublefastaref_out}", mode: 'copy'
  }
  
  input:
    each(fasta)

  output:
    tuple val(params.hbvdb), path("doubled/${fasta.baseName}_doubled.fasta"), emit: doubledfasta

  script:
    """
    mkdir doubled
    seqkit concat ${fasta} ${fasta} -o doubled/${fasta.baseName}_doubled.fasta
    """
}
