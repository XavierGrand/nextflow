version = "407"
container_url = "lbmc/ucsc:${version}"

include {
  index_fasta
} from './../samtools/main'

params.bedgraph_to_bigwig = ""
params.bedgraph_to_bigwig_out = ""
process bedgraph_to_bigwig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bedgraph_to_bigwig_out != "") {
    publishDir "results/${params.bedgraph_to_bigwig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bg)
  tuple val(file_id), path(bed)

  output:
  tuple val(file_id), path("*.bw"), emit: bw

  script:
"""
LC_COLLATE=C
# transform bed file of start-stop chromosome size to stop chromosome size
awk -v OFS="\\t" '{print \$1, \$3}' ${bed} > chromsize.txt

sort -T ./ -k1,1 -k2,2n ${bg} > \
  bedGraphToBigWig ${params.bedgraph_to_bigwig} - \
    chromsize.txt \
    ${bg.simpleName}_norm.bw
"""
}

params.wig_to_bedgraph = ""
params.wig_to_bedgraph_out = ""
workflow wig_to_bedgraph {
  take:
    fasta
    wig
  main:
    wig_to_bigwig(
      fasta,
      wig
    )
    bigwig_to_bedgraph(
      wig_to_bigwig.out.bw
    )
  emit:
  bg = bigwig_to_bedgraph.out.bg
}

workflow wig2_to_bedgraph2 {
  take:
    fasta
    wig
  main:
    wig2_to_bigwig2(
      fasta,
      wig
    )
    bigwig2_to_bedgraph2(
      wig2_to_bigwig2.out.bw
    )
  emit:
  bg = bigwig2_to_bedgraph2.out.bg
}

params.bigwig_to_bedgraph = ""
params.bigwig_to_bedgraph_out = ""
process bigwig_to_bedgraph {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_bedgraph_out != "") {
    publishDir "results/${params.bigwig_to_bedgraph_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bw)

  output:
  tuple val(file_id), path("*.bg"), emit: bg

  script:
"""
bigWigToBedGraph ${bw} ${bw.simpleName}.bg
"""
}

params.bigwig2_to_bedgraph2 = ""
params.bigwig2_to_bedgraph2_out = ""
process bigwig2_to_bedgraph2 {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_bedgraph_out != "") {
    publishDir "results/${params.bigwig_to_bedgraph_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bw_a), path(bw_b)

  output:
  tuple val(file_id), path("${bw_a.simpleName}.bg"), path("${bw_b.simpleName}.bg"), emit: bg

  script:
"""
bigWigToBedGraph ${bw_a} ${bw_a.simpleName}.bg
bigWigToBedGraph ${bw_b} ${bw_b.simpleName}.bg
"""
}

params.bigwig_to_wig = ""
params.bigwig_to_wig_out = ""
process bigwig_to_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bw)

  output:
  tuple val(file_id), path("*.wig"), emit: wig

  script:
"""
bigWigToBedGraph ${bw} ${bw.simpleName}.bg
bedgraph_to_wig.pl --bedgraph ${bw.simpleName}.bg --wig ${bw.simpleName}.wig --step 10
"""
}

params.bigwig2_to_wig2 = ""
params.bigwig2_to_wig2_out = ""
process bigwig2_to_wig2 {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(bw_a), path(bw_b)

  output:
  tuple val(file_id), path("${bw_a.simpleName}.wig"), path("${bw_b.simpleName}.wig"), emit: wig

  script:
"""
bigWigToBedGraph ${bw_a} ${bw_a.simpleName}.bg
bedgraph_to_wig.pl --bedgraph ${bw_a.simpleName}.bg --wig ${bw_a.simpleName}.wig --step 10
bigWigToBedGraph ${bw_b} ${bw_b.simpleName}.bg
bedgraph_to_wig.pl --bedgraph ${bw_b.simpleName}.bg --wig ${bw_b.simpleName}.wig --step 10
"""
}

params.wig_to_bigwig = ""
params.wig_to_bigwig_out = ""

workflow wig_to_bigwig {
  take:
    fasta
    wig
  main:
    index_fasta(fasta)
    wig_to_bigwig_sub(
      wig,
      index_fasta.out.index
    )
  emit:
  bw = wig_to_bigwig_sub.out.bw
}

process wig_to_bigwig_sub {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(w)
  tuple val(idx_id), path(fasta_idx)

  output:
  tuple val(file_id), path("${w.simpleName}.bw"), emit: bw

  script:
"""
cut -f 1,2 ${fasta_idx} > ${fasta_idx.simpleName}.sizes
wigToBigWig -clip ${w} ${fasta_idx.simpleName}.sizes ${w.simpleName}.bw
"""
}

params.wig2_to_bigwig2 = ""
params.wig2_to_bigwig2_out = ""

workflow wig2_to_bigwig2 {
  take:
    fasta
    wigs
  main:
    index_fasta(fasta)
    wig2_to_bigwig2_sub(
      wigs,
      index_fasta.out.index
    )
  emit:
  bw = wig2_to_bigwig2_sub.out.bw
}

process wig2_to_bigwig2_sub {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "${file_id}"
  if (params.bigwig_to_wig_out != "") {
    publishDir "results/${params.bigwig_to_wig_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(w_a), path(w_b)
  tuple val(idx_id), path(fasta_idx)

  output:
  tuple val(file_id), path("${w_a.simpleName}.bw"), path("${w_b.simpleName}.bw"), emit: bw

  script:
"""
cut -f 1,2 ${fasta_idx} > ${fasta_idx.simpleName}.sizes
wigToBigWig -clip ${w_a} ${fasta_idx.simpleName}.sizes ${w_a.simpleName}.bw
wigToBigWig -clip ${w_b} ${fasta_idx.simpleName}.sizes ${w_b.simpleName}.bw
"""
}