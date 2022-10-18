version = "1.11"
container_url = "lbmc/samtools:${version}"

params.index_fasta = ""
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
  output:
    tuple val(file_id), path("*.fai"), emit: index

  script:
"""
if gzip -t ${fasta}; then
  zcat ${fasta} > ${fasta.simpleName}.fasta
  samtools faidx ${params.index_fasta}  ${fasta.simpleName}.fasta
else
  samtools faidx ${params.index_fasta} ${fasta}
fi

"""
}

params.filter_bam_quality_threshold = 30
params.filter_bam_quality = "-q ${params.filter_bam_quality_threshold}"
params.filter_bam_quality_out = ""
process filter_bam_quality {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.filter_bam_quality_out != "") {
    publishDir "results/${params.filter_bam_quality_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_filtered.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} ${params.filter_bam_quality} > \
  ${bam.simpleName}_filtered.bam
"""
}

params.filter_bam = ""
params.filter_bam_out = ""
process filter_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.filter_bam_out != "") {
    publishDir "results/${params.filter_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)
    tuple val(bed_id), path(bed)

  output:
    tuple val(file_id), path("*_filtered.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} -hb ${bam} -L ${bed} ${params.filter_bam} > \
  ${bam.simpleName}_filtered.bam
"""
}

params.rm_from_bam = ""
params.rm_from_bam_out = ""
process rm_from_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.rm_from_bam_out != "") {
    publishDir "results/${params.rm_from_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)
    tuple val(bed_id), path(bed)

  output:
    tuple val(file_id), path("*_filtered.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} ${params.filter_bam} -hb -L ${bed} -U ${bam.simpleName}_filtered.bam ${bam} >  /dev/null
"""
}

params.filter_bam_mapped = "-F 4"
params.filter_bam_mapped_out = ""
process filter_bam_mapped {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.filter_bam_mapped_out != "") {
    publishDir "results/${params.filter_bam_mapped_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_mapped.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} ${params.filter_bam_mapped} -hb ${bam} > \
  ${bam.simpleName}_mapped.bam
"""
}

params.filter_bam_unmapped = "-f 4"
params.filter_bam_unmapped_out = ""
process filter_bam_unmapped {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.filter_bam_unmapped_out != "") {
    publishDir "results/${params.filter_bam_unmapped_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_unmapped.bam"), emit: bam
  script:
"""
samtools view -@ ${task.cpus} ${params.filter_bam_unmapped} -hb ${bam} > ${bam.simpleName}_unmapped.bam
"""
}

params.index_bam = ""
params.index_bam_out = ""
process index_bam {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_bam_out != "") {
    publishDir "results/${params.index_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("${bam}"), path("*.bam.bai"), emit: bam_idx

  script:
"""
samtools index ${params.index_bam} ${bam}
"""
}

params.sort_bam = ""
params.sort_bam_out = ""
process sort_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.sort_bam_out != "") {
    publishDir "results/${params.sort_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
samtools sort -@ ${task.cpus} ${params.sort_bam} -O BAM -o ${bam.simpleName}_sorted.bam ${bam}
"""
}

params.split_bam = ""
params.split_bam_out = ""
params.flag_value = "0x10"
process split_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.split_bam_out != "") {
    publishDir "results/${params.split_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_forward.bam*"), emit: bam_forward
    tuple val(file_id), path("*_reverse.bam*"), emit: bam_reverse
  script:
"""
samtools view -@ ${Math.round(task.cpus/2)} ${params.split_bam} \
  -hb -F ${flag_value} ${bam} > ${bam.simpleName}_forward.bam &
samtools view -@ ${Math.round(task.cpus/2)} ${params.split_bam} \
  -hb -f ${flag_value} ${bam} > ${bam.simpleName}_reverse.bam
"""
}

params.merge_bam = ""
params.merge_bam_out = ""
process merge_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.merge_bam_out != "") {
    publishDir "results/${params.merge_bam_out}", mode: 'copy'
  }

  input:
    tuple val(first_file_id), path(first_bam)
    tuple val(second_file_id), path(second_bam)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam
  script:
"""
samtools merge -@ ${task.cpus} ${params.merge_bam} ${first_bam} ${second_bam} \
  ${first_bam.simpleName}_${second_file.simpleName}.bam
"""
}

params.merge_multi_bam = ""
params.merge_multi_bam_out = ""
process merge_multi_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.merge_multi_bam_out != "") {
    publishDir "results/${params.merge_multi_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bams)

  output:
    tuple val(file_id), path("*_merged.bam*"), emit: bam
  script:
"""
samtools merge -@ ${task.cpus} \
  ${params.merge_multi_bam} \
  ${bams[0].simpleName}_merged.bam \
  ${bams}
"""
}

params.stats_bam = ""
params.stats_bam_out = ""
process stats_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.stats_bam_out != "") {
    publishDir "results/${params.stats_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*.tsv"), emit: tsv
    path "*.flagstat.txt", emit: report 
  script:
"""
samtools flagstat -@ ${task.cpus} ${params.stats_bam} -O tsv ${bam} > ${bam.simpleName}.flagstat.txt
cp ${bam.simpleName}.flagstat.txt ${bam.simpleName}.tsv
"""
}

params.flagstat_2_multiqc = ""
params.flagstat_2_multiqc_out = ""
process flagstat_2_multiqc {
  tag "$file_id"
  if (params.flagstat_2_multiqc_out != "") {
    publishDir "results/${params.flagstat_2_multiqc_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(tsv)

  output:
    tuple val(file_id), path("*.txt"), emit: report
"""
mv ${tsv} ${tsv.simpleName}.flagstat.txt
"""
}

params.idxstat_2_multiqc = ""
params.idxstat_2_multiqc_out = ""
process idxstat_2_multiqc {
  tag "$file_id"
  if (params.idxstat_2_multiqc_out != "") {
    publishDir "results/${params.idxstat_2_multiqc_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(tsv)

  output:
    tuple val(file_id), path("*.txt"), emit: report
"""
mv ${tsv} ${tsv.simpleName}.idxstats.txt
"""
}

params.split_rt_bam = ""
params.split_rt_bam_out = ""
params.flag_fwd = "99,147"
params.flag_rev = "83,163"
process split_rt_bam {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.split_rt_bam_out != "") {
    publishDir "results/${params.split_rt_bam_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)

  output:
    tuple val(file_id), path("*_forward.bam*"), emit: bam_forward
    tuple val(file_id), path("*_reverse.bam*"), emit: bam_reverse

  script:
"""
samtools view -@ ${Math.round(task.cpus/2)} ${params.split_rt_bam} \
  -hb -F${params.flag_fwd} ${bam} > ${bam.simpleName}_forward.bam &
samtools view -@ ${Math.round(task.cpus/2)} ${params.split_rt_bam} \
  -hb -F${params.flag_rev} ${bam} > ${bam.simpleName}_reverse.bam
"""
}