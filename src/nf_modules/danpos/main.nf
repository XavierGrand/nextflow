version = "v2.2.2_cv3"
container_url = "biocontainers/danpos:${version}"

include {
  bigwig2_to_wig2;
  bigwig_to_wig;
  wig_to_bedgraph;
  wig2_to_bedgraph2
} from "./../ucsc/main.nf"

params.dpos = "--smooth_width 0 -n N "
params.dpos_out = ""

workflow dpos_bam_bg {
  take:
    fasta
    fastq
    bam

  main:
    dpos_bam(fastq, bam)
    wig2_to_bedgraph2(fasta, dpos_bam.out.wig)

  emit:
    bg = wig2_to_bedgraph2.out.bg
    wig = dpos_bam.out.wig
    bed = dpos_bam.out.bed
}

process dpos_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id), path(bam_ip), path(bam_wce)

  output:
    tuple val(file_id), path("${file_prefix}/${bam_ip.simpleName}*.wig"), path("${file_prefix}/${bam_wce.simpleName}*.wig"), emit: wig
  tuple val(file_id), path("${file_prefix}/*.positions.bed"), emit: bed

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -b ${bam_wce} \
  -o ${file_prefix} \
  ${bam_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.bed
"""
}

workflow dpos_bw {
  take:
    fasta
    fastq
    bw
  main:
    bigwig2_to_wig2(bw)
    dpos_wig(fastq, bigwig2_to_wig2.out.wig)
    wig_to_bedgraph(fasta, bigwig2_to_wig2.out.wig)

  emit:
  bg = wig_to_bedgraph.out.bg
  wig = bigwig2_to_wig2.out.wig
  bed = dpos_wig.out.bed
}

process dpos_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id), path(wig_ip), path(wig_wce)

  output:
    tuple val(file_id), path("${file_prefix}/*.positions.bed"), emit: bed
    tuple val(file_id), path("${file_prefix}/${bam_ip.simpleName}*.wig"), path("${file_prefix}/${bam_wce.simpleName}*.wig"), emit: wig

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -b ${wig_wce} \
  -o ${file_prefix} \
  ${wig_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${wig_ip.simpleName}.positions.xls > ${file_prefix}/${wig_ip.simpleName}.positions.bed
"""
}

workflow dpos_bw_no_b {
  take:
    fasta
    fastq
    bw
  main:
    bigwig_to_wig(bw)
    dpos_wig_no_b(fastq, bigwig_to_wig.out.wig)
    wig_to_bedgraph(fasta, bigwig_to_wig.out.wig)

  emit:
  bg = wig_to_bedgraph.out.bg
  wig = bigwig_to_wig.out.wig
  bed = dpos_wig_no_b.out.bed
}

process dpos_wig_no_b {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id), path(wig_ip)

  output:
    tuple val(file_id), path("${file_prefix}/*.positions.bed"), emit: bed

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -o ${file_prefix} \
  ${wig_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${wig_ip.simpleName}.positions.xls > ${file_prefix}/${wig_ip.simpleName}.positions.bed
"""
}

workflow dwig_bwvsbw {
  take:
    fasta
    fastq
    bw_a
    bw_b
  main:
    dpos_wigvswig(
      fastq,
      bigwig2_to_wig2(bw_a),
      bigwig2_to_wig2(bw_b),
    )
    wig_to_bedgraph(fasta, dpos_wigvswig.out.wig)

  emit:
  bg = wig_to_bedgraph.out.bg
  wig = dpeak_wig.out.wig
  bed = dpeak_wig.out.bed
}

process dpos_wigvswig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpos_out != "") {
    publishDir "results/${params.dpos_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id_a), path(wig_ip_a)
    tuple val(file_id_b), path(wig_ip_b)

  output:
    tuple val(file_id), path("${file_prefix}/${wig_ip_a.simpleName}*.wig"), emit: wig
  tuple val(file_id), path("${file_prefix}/*.positions.bed"), emit: bed

  script:

  switch(file_id_a) {
    case {it instanceof List}:
      file_prefix = file_id_a[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id_a.values()[0]
    break
    default:
      file_prefix = file_id_a
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpos -m ${m} \
  ${params.dpos} \
  -b ${wig_ip_a},${wig_ip_b} \
  -o ${file_prefix} \
  ${wig_ip_a}:${wig_ip_b}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.positions.xls > ${file_prefix}/${bam_ip.simpleName}.positions.bed
"""
}

params.dpeak = "--smooth_width 0 -n N "
params.dpeak_out = ""

process dpeak_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id), path(bam_ip), path(bam_wce)

  output:
    tuple val(file_id), path("${file_prefix}/${bam_ip.simpleName}*.wig"), path("${file_prefix}/${bam_wce.simpleName}*.wig"), emit: wig
  tuple val(file_id), path("${file_prefix}/*.positions.bed"), path("${file_prefix}/*.summit.bed"), emit: bed
    tuple val(file_id), path("${file_prefix}/*.bed"), emit: bed

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpeak -m ${m} \
  ${params.dpeak} \
  -b ${bam_wce} \
  -o ${file_prefix} \
  ${bam_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.bed
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$4-1, \$4, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.summit.bed
"""
}

workflow dpeak_bw {
  take:
    fasta
    fastq
    bw
  main:
    dpeak_wig(fastq, bigwig2_to_wig2(bw))
    wig2_to_bedgraph2(fasta, dpeak_wig.out.wig)

  emit:
  bg = wig2_to_bedgraph2.out.bg
  wig = dpeak_wig.out.wig
  bed = dpeak_wig.out.bed
}


process dpeak_wig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id), path(wig_ip), path(wig_wce)

  output:
  tuple val(file_id), path("${file_prefix}/${wig_ip.simpleName}.bgsub.wig"), path("${file_prefix}/${wig_wce.simpleName}.wig"), emit: wig
  tuple val(file_id), path("${file_prefix}/*.positions.bed"), path("${file_prefix}/*.summit.bed"), emit: bed

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpeak -m ${m} \
  ${params.dpeak} \
  -b ${wig_wce} \
  -o ${file_prefix} \
  ${wig_ip}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${wig_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${wig_ip.simpleName}.bgsub.positions.bed
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$4-1, \$4, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${wig_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${wig_ip.simpleName}.bgsub.positions.summit.bed
"""
}

workflow dpeak_bwvsbw {
  take:
    fasta
    fastq
    bw_a
    bw_b
  main:
    dpeak_wigvswig(
      fastq,
      bigwig2_to_wig2(bw_a),
      bigwig2_to_wig2(bw_b),
    )
    wig2_to_bedgraph2(fasta, dpeak_wigvswig.out.wig)

  emit:
  bg = wig2_to_bedgraph2.out.bg
  wig = dpeak_wig.out.wig
  bed = dpeak_wig.out.bed
}


process dpeak_wigvswig {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.dpeak_out != "") {
    publishDir "results/${params.dpeak_out}", mode: 'copy', overwrite: true
  }

  input:
    val fastq 
    tuple val(file_id_a), path(wig_ip_a), path(wig_wce_a)
    tuple val(file_id_b), path(wig_ip_b), path(wig_wce_b)

  output:
  tuple val(file_id), path("${file_prefix}/${wig_ip_a.simpleName}.bgsub.wig"), path("${file_prefix}/${wig_wce_a.simpleName}.wig"), emit: wig
  tuple val(file_id), path("${file_prefix}/*.positions.bed"), path("${file_prefix}/*.summit.bed"), emit: bed

  script:

  switch(file_id_a) {
    case {it instanceof List}:
      file_prefix = file_id_a[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id_a.values()[0]
    break
    default:
      file_prefix = file_id_a
    break
  }

  m = 0
  if (fastq[1].size() == 2){
    m = 1
  }
"""
danpos.py dpeak -m ${m} \
  ${params.dpeak} \
  -b ${wig_ip_a}:${wig_wce_a},${wig_ip_b}:${wig_wce_b} \
  -o ${file_prefix} \
  ${wig_ip_a}:${wig_ip_b}
mv ${file_prefix}/pooled/* ${file_prefix}/
rm -R ${file_prefix}/pooled
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$2-1, \$3, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.bed
awk -v FS='\t' -v OFS='\t' 'FNR > 1 { print \$1, \$4-1, \$4, "Interval_"NR-1, \$6, "+" }' ${file_prefix}/${bam_ip.simpleName}.bgsub.peaks.xls > ${file_prefix}/${bam_ip.simpleName}.bgsub.positions.summit.bed
"""
}