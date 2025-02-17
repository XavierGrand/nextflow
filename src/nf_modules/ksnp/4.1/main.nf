version = "4.1"
container_url = "xgrand/ksnp:${version}"

params.mkinfile = ""
params.mkinfile_out = ""
process mkinfile {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "$file_id"
  if (params.mkinfile_out != "") {
    publishDir "results/${params.mkinfile_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(indir)
  output:
    tuple val(file_id), path("*.input"), emit: infile

  script:
"""
MakeKSNP4infile -indir ${indir} -outfile KSNP4.input
"""
}

params.mktree = "-k 11 -vcf -ML"
params.mktree_out = ""
process mktree {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.mktree_out != "") {
    publishDir "results/${params.mktree_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(infile)

  script:
"""
kSNP4 ${params.mktree} -CPU ${task.cpus} -in ${infile} -outdir . | tee kSNP4.log
"""
}

params.kmersize = ""
params.kmersize_out = ""
process kmersize {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.kmersize_out != "") {
    publishDir "results/${params.kmersize_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(infile)

  script:
"""
Kchooser4 ${params.kmersize} -in ${infile} > ideal_kmer_size.txt
"""
}