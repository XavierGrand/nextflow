version = "4.2.0.0"
container_url = "broadinstitute/gatk:${version}"

def get_file_prefix(file_id) {
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else if (file_id instanceof Map) {
      library = file_id[0]
      file_prefix = file_id[0]
      if (file_id.containsKey('library')) {
        library = file_id.library
        file_prefix = file_id.id
      }
  } else {
    file_prefix = file_id
  }
  return file_prefix
}

include {
  index_fasta as samtools_index_fasta;
  index_bam;
} from './../samtools/main.nf'
include {
  index_fasta as picard_index_fasta;
  index_bam as picard_index_bam;
  mark_duplicate;
} from './../picard/main.nf'

params.variant_calling_out = ""
workflow germline_cohort_data_variant_calling {
  take:
    bam
    fasta
  main:
    // data preparation
    mark_duplicate(bam)
    index_bam(mark_duplicate.out.bam)
    picard_index_bam(mark_duplicate.out.bam)
    index_bam.out.bam_idx
      .join(picard_index_bam.out.index)
      .set{ bam_idx }
    picard_index_fasta(fasta)
    samtools_index_fasta(fasta)
    fasta
      .join(picard_index_fasta.out.index)
      .join(samtools_index_fasta.out.index)
      .set{ fasta_idx }
    
    // variant calling
    call_variants_per_sample(
      bam_idx,
      fasta_idx.collect()
    )
    call_variants_all_sample(
      call_variants_per_sample.out.gvcf,
      fasta_idx
    )
  emit:
    vcf = call_variants_all_sample.out.vcf
}

/*******************************************************************/
workflow base_quality_recalibrator{
  take:
    bam_idx
    fasta_idx
    vcf

  main:
    index_vcf(vcf)
    compute_base_recalibration(
      bam_idx,
      fasta_idx,
      index_vcf.out.vcf_idx
    ) 
    apply_base_recalibration(
      bam_idx,
      fasta_idx,
      compute_base_recalibration.out.table
    )
    emit:
    bam = apply_base_recalibration.out.bam
}

process index_vcf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  input:
    tuple val(file_id), path(vcf)
  output:
    tuple val(file_id), path("${vcf}"), path("*"), emit: vcf_idx

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${vcf}
"""
}

process compute_base_recalibration {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  input:
    tuple val(file_id), path(bam), path(bam_idx), path(bam_idx_bis)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
    tuple val(vcf_id), path(vcf), path(vcf_idx)
  output:
    tuple val(file_id), path("${bam.simpleName}.table"), emit: table

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
  def vcf_cmd = ""
  if (vcf instanceof List){
    for (vcf_file in vcf){
      vcf_cmd += "--known-sites ${vcf_file} "
    }
  } else {
    vcf_cmd = "--known-sites ${vcf} "
  }
"""
 gatk --java-options "-Xmx${xmx_memory}G" BaseRecalibrator \
   -I ${bam} \
   -R ${fasta} \
   ${vcf_cmd} \
   -O ${bam.simpleName}.table
"""
}

process apply_base_recalibration {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  input:
    tuple val(file_id), path(bam), path(bam_idx), path(bam_idx_bis)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
    tuple val(table_id), path(table)
  output:
    tuple val(file_id), path("${bam.simpleName}_recalibrate.bam"), emit: bam

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
 gatk --java-options "-Xmx${xmx_memory}G" ApplyBQSR \
   -R ${fasta} \
   -I ${bam} \
   --bqsr-recal-file ${table} \
   -O ${bam.simpleName}_recalibrate.bam
"""
}

/*******************************************************************/
params.variant_calling_gvcf_out = ""
process call_variants_per_sample {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.variant_calling_gvcf_out != "") {
    publishDir "results/${params.variant_calling_gvcf_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(bam), path(bam_idx), path(bam_idx_bis)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("${bam.simpleName}.gvcf.gz"), emit: gvcf

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
 gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller  \
   -R ${fasta} \
   -I ${bam} \
   -O ${bam.simpleName}.gvcf.gz \
   -ERC GVCF
"""
}

/*******************************************************************/

workflow call_variants_all_sample {
  take:
    gvcf
    fasta_idx

  main:
    index_gvcf(gvcf)
    validate_gvcf(
      index_gvcf.out.gvcf_idx,
      fasta_idx.collect()
    )
    consolidate_gvcf(
      validate_gvcf.out.gvcf
      .groupTuple(),
      fasta_idx.collect()
    )
    genomic_db_call(
      consolidate_gvcf.out.gvcf_idx,
      fasta_idx.collect()
    )
  emit:
    vcf = genomic_db_call.out.vcf
}

process index_gvcf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  input:
    tuple val(file_id), path(gvcf)
  output:
    tuple val(file_id), path("${gvcf}"), path("${gvcf}.tbi"), emit: gvcf_idx
    tuple val(file_id), path("${gvcf.simpleName}_IndexFeatureFile_report.txt"), emit: report

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
      -I ${gvcf} 2> ${gvcf.simpleName}_IndexFeatureFile_report.txt
"""
}

process validate_gvcf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  input:
    tuple val(file_id), path(gvcf), path(gvcf_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("${gvcf}"), path("${gvcf_idx}"), emit: gvcf

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" ValidateVariants \
   -V ${gvcf} \
   -R ${fasta} -gvcf
"""
}

process consolidate_gvcf {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  input:
    tuple val(file_id), path(gvcf), path(gvcf_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("${file_prefix}.gvcf"), path("${file_prefix}.gvcf.idx"), emit: gvcf_idx
    tuple val(file_id), path("${file_prefix}_CombineGVCFs_report.txt"), emit: report

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
  def gvcf_cmd = ""
  if (gvcf instanceof List){
    for (gvcf_file in gvcf){
      gvcf_cmd += "-V ${gvcf_file} "
    }
  } else {
    gvcf_cmd = "-V ${gvcf} "
  }
"""
mkdir tmp
gatk --java-options "-Xmx${xmx_memory}G" CombineGVCFs \
    ${gvcf_cmd} \
    -R ${fasta} \
    -O ${file_prefix}.gvcf 2> ${file_prefix}_CombineGVCFs_report.txt
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
      -I ${file_prefix}.gvcf 2> ${file_prefix}_IndexFeatureFile_report.txt
"""
}

process genomic_db_call {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.variant_calling_out != "") {
    publishDir "results/${params.variant_calling_out}", mode: 'copy'
  }
  input:
    tuple val(file_id), path(gvcf), path(gvcf_idx)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("${gvcf.simpleName}.vcf.gz"), emit: vcf

  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
  def gvcf_cmd = ""
  if (gvcf instanceof List){
    for (gvcf_file in gvcf){
      gvcf_cmd += "--V ${gvcf_file} "
    }
  } else {
    gvcf_cmd = "--V ${gvcf} "
  }
"""
mkdir tmp
gatk --java-options "-Xmx${xmx_memory}G" GenotypeGVCFs \
   -R ${fasta} \
   -V ${gvcf} \
   -O ${gvcf.simpleName}.vcf.gz \
   --tmp-dir ./tmp
"""
}

/*******************************************************************/
params.variant_calling = ""
process variant_calling {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller \
  ${params.variant_calling} \
  -R ${fasta} \
  -I ${bam} \
  -O ${bam.simpleName}.vcf
"""
}

params.filter_snp = ""
params.filter_snp_out = ""
process filter_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  ${params.filter_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type SNP \
  -O ${vcf.simpleName}_snp.vcf
"""
}

params.filter_indels = ""
params.filter_indels_out = ""
process filter_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  ${params.filter_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type INDEL \
  -O ${vcf.simpleName}_indel.vcf
"""
}

params.high_confidence_snp_filter = "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)"
params.high_confidence_snp = "--filter-expression \"${params.high_confidence_snp_filter}\" --filter-name \"basic_snp_filter\""
params.high_confidence_snp_out = ""
process high_confidence_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" VariantFiltration \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_snp} \
  -O ${vcf.simpleName}_filtered_snp.vcf
"""
}

params.high_confidence_indel_filter = "QD < 3.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0"
params.high_confidence_indels = "--filter-expression \"${params.high_confidence_indel_filter}\" --filter-name \"basic_indel_filter\""
params.high_confidence_indels_out = ""
process high_confidence_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" VariantFiltration \
  -R ${fasta} \
  -V ${vcf} \
  ${params.high_confidence_indels} \
  -O ${vcf.simpleName}_filtered_indel.vcf
"""
}

params.recalibrate_snp_table = ""
params.recalibrate_snp_table_out = ""
process recalibrate_snp_table {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.recalibrate_snp_table_out != "") {
    publishDir "results/${params.recalibrate_snp_table_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx), path(bam_idx_bis)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("recal_data_table"), emit: recal_table
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${snp_file}
gatk --java-options "-Xmx${xmx_memory}G" IndexFeatureFile \
  -I ${indel_file}
gatk --java-options "-Xmx${xmx_memory}G" BaseRecalibrator \
  ${params.recalibrate_snp_table} \
  -R ${fasta} \
  -I ${bam} \
  -known-sites ${snp_file} \
  -known-sites ${indel_file} \
  -O recal_data_table
"""
}

params.recalibrate_snp = ""
params.recalibrate_snp_out = ""
process recalibrate_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.recalibrate_snp_out != "") {
    publishDir "results/${params.recalibrate_snp_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(snp_file), path(indel_file), path(bam), path(bam_idx), path(recal_table)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.bam"), emit: bam
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" ApplyBQSR \
  ${params.recalibrate_snp} \
  -R ${fasta} \
  -I ${bam} \
  --bqsr-recal-file recal_data_table \
  -O ${bam.simpleName}_recal.bam
"""
}

params.haplotype_caller = ""
params.haplotype_caller_out = ""
process haplotype_caller {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" HaplotypeCaller \
  ${params.haplotype_caller} \
  -R ${fasta} \
  -I ${bam} \
  -ERC GVCF \
  -O ${bam.simpleName}.gvcf
"""
}

params.gvcf_genotyping = ""
params.gvcf_genotyping_out = ""
process gvcf_genotyping {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.gvcf_genotyping_out != "") {
    publishDir "results/${params.gvcf_genotyping_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(gvcf)
    tuple val(ref_id), path(fasta), path(fai), path(dict)
  output:
    tuple val(file_id), path("*.vcf.gz"), emit: vcf
  script:
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" GenotypeGVCFs \
  ${params.gvcf_genotyping} \
  -R ${fasta} \
  -V ${gvcf} \
  -O ${gvcf.simpleName}_joint.vcf.gz
"""
}

params.select_variants_snp = ""
params.select_variants_snp_out = ""
process select_variants_snp {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}GG" SelectVariants \
  ${params.select_variants_snp} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type SNP \
  -O ${vcf.simpleName}_joint_snp.vcf
"""
}

params.select_variants_indels = ""
params.select_variants_indels_out = ""
process select_variants_indels {
  container = "${container_url}"
  label "big_mem_mono_cpus"
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" SelectVariants \
  ${params.select_variants_indels} \
  -R ${fasta} \
  -V ${vcf} \
  -select-type INDEL \
  -O ${file_prefix}_joint_indel.vcf
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
  xmx_memory = "${task.memory}" - ~/\s*GB/
  file_prefix = get_file_prefix(file_id)
"""
gatk --java-options "-Xmx${xmx_memory}G" FastaAlternateReferenceMaker\
  ${params.personalized_genome} \
  -R ${reference} \
  -V ${vcf} \
  -O ${vcf.simpleName}_genome.fasta
"""
}



