version = "5.1_24Aug19.3e8--hdfd78af_1"
container_url = "quay.io/biocontainers/beagle::${version}"

params.phasing = ""
process phasing {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"

  input:
    tuple val(file_id), path(vcf)
    tuple val(ref_id), path(ref_vcf)

  output:
    tuple val(file_id), path("*.bam*"), emit: bam

  script:
"""
beagle nthread=${task.cpus} \
  gtgl=${vcf} \
  ref=${ref_vcf}
"""
}
