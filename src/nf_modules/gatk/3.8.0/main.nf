version = "3.8.0"
container_url = "lbmc/gatk:${version}"

params.variant_calling = ""
params.variant_calling_out = ""
process variant_calling {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.variant_calling_out != "") {
    publishDir "results/${params.variant_calling_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam), path(bai)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T HaplotypeCaller \
  -nct ${task.cpus} \
  ${params.variant_calling} \
  -R ${fasta} \
  -I ${bam} \
  -o ${file_prefix}.vcf
"""
}

params.filter_snp = ""
params.filter_snp_out = ""
process filter_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.filter_snp_out != "") {
    publishDir "results/${params.filter_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.filter_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType SNP \
  -o ${file_prefix}_snp.vcf
"""
}

params.filter_indels = ""
params.filter_indels_out = ""
process filter_indels {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.filter_indels_out != "") {
    publishDir "results/${params.filter_indels_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.filter_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType INDEL \
  -o ${file_prefix}_indel.vcf
"""
}

params.high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"
params.high_confidence_snp = "--filterExpression \"${params.high_confidence_snp_filter}\" --filterName \"basic_snp_filter\""
params.high_confidence_snp_out = ""
process high_confidence_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.high_confidence_snp_out != "") {
    publishDir "results/${params.high_confidence_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_snp.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_snp} \
  -o ${file_prefix}_filtered_snp.vcf
"""
}

params.high_confidence_indel_filter = "QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"
params.high_confidence_indels = "--filterExpression \"${params.high_confidence_indel_filter}\" --filterName \"basic_indel_filter\""
params.high_confidence_indels_out = ""
process high_confidence_indels {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.high_confidence_indels_out != "") {
    publishDir "results/${params.high_confidence_indels_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_indel.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T VariantFiltration \
  -nct ${task.cpus} \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_indels} \
  -o ${file_prefix}_filtered_indel.vcf
"""
}

params.recalibrate_snp_table = ""
params.recalibrate_snp_table_out = ""
process recalibrate_snp_table {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.recalibrate_snp_table_out != "") {
    publishDir "results/${params.recalibrate_snp_table_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("recal_data_table"), emit: recal_table
  script:
"""
gatk3 -T BaseRecalibrator \
  -nct ${task.cpus} \
  ${recalibrate_snp_table} \
  -R ${fasta} \
  -I ${bam} \
  -knownSites ${snp_file} \
  -knownSites ${indel_file} \
  -o recal_data_table
"""
}

params.recalibrate_snp = ""
params.recalibrate_snp_out = ""
process recalibrate_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.recalibrate_snp_out != "") {
    publishDir "results/${params.recalibrate_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx)
    tuple val(table_id), path(recal_data_table)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.bam"), emit: bam
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T PrintReads \
  --use_jdk_deflater \
  --use_jdk_inflater \
  ${recalibrate_snp} \
  -nct ${task.cpus} \
  -R ${fasta} \
  -I ${bam} \
  -BQSR recal_data_table \
  -o ${file_prefix}_recal.bam
"""
}

params.haplotype_caller = ""
params.haplotype_caller_out = ""
process haplotype_caller {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.haplotype_caller_out != "") {
    publishDir "results/${params.haplotype_caller_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.gvcf"), emit: gvcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T HaplotypeCaller \
  -nct ${task.cpus} \
  ${params.haplotype_caller} \
  -R ${fasta} \
  -I ${bam} \
  -ERC GVCF \
  -variant_index_type LINEAR -variant_index_parameter 128000 \
  -o ${file_prefix}.gvcf
"""
}

params.gvcf_genotyping = ""
params.gvcf_genotyping_out = ""
process gvcf_genotyping {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.gvcf_genotyping_out != "") {
    publishDir "results/${params.gvcf_genotyping_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gvcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T GenotypeGVCFs \
  -nct ${task.cpus} \
  ${params.gvcf_genotyping} \
  -R ${fasta} \
  -V ${gvcf} \
  -o ${file_prefix}_joint.vcf
"""
}

params.select_variants_snp = ""
params.select_variants_snp_out = ""
process select_variants_snp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.select_variants_snp_out != "") {
    publishDir "results/${params.select_variants_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_snp.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.select_variants_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType SNP \
  -o ${file_prefix}_joint_snp.vcf
"""
}

params.select_variants_indels = ""
params.select_variants_indels_out = ""
process select_variants_indels {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.select_variants_indels_out != "") {
    publishDir "results/${params.select_variants_indels_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_joint_indel.vcf"), emit: vcf
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T SelectVariants \
  -nct ${task.cpus} \
  ${params.select_variants_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -selectType INDEL \
  -o ${file_prefix}_joint_indel.vcf
"""
}

params.personalized_genome = ""
params.personalized_genome_out = ""
process personalized_genome {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.personalized_genome_out != "") {
    publishDir "results/${params.personalized_genome_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*_genome.fasta"), emit: fasta

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
"""
gatk3 -T FastaAlternateReferenceMaker\
  ${params.personalized_genome} \
  -R ${reference} \
  -V ${vcf} \
  -o ${file_prefix}_genome.fasta
"""
}

