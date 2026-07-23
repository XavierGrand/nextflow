version = "0.8.0"
container_url = "xgrand/dorado:${version}"

params.basecalling_out = ""
params.bc_kit = "" // Barcoding kit : "EXP-PBC001"
params.min_qscore = 8
params.cuda = "cuda:all"
params.models = "sup"

process basecall_fast5 {
  container = "${container_url}"
  label "gpus"
  tag "$fast5_folder"
  if (params.basecalling_out != "") {
    publishDir "results/${params.basecalling_out}", mode: 'copy'
  }

  input:
    path(pod5)

  output:
    path("*.fastq"), emit: fastq
    path("barcoding_summary.txt"), emit: barcoding_summary

  script:
  // Cas 1: Input fast5 BARCODED
  if (params.bc_kit != "")
"""
dorado basecaller \
	-x ${params.cuda} \
	--recursive \
	--min-qscore ${params.min_qscore} \
	--output-dir basecalled/ \
  --kit-name ${params.bc_kit} \
  --no-trim \
	${params.models} \
	${pod5}

dorado demux \
  --emit-fastq \
  --emit-summary \
  --threads ${task.cpus} \
  --kit-name ${params.bc_kit} \
  --output-dir . \
  basecalled/
"""
  else
"""
dorado basecaller \
	-x ${params.cuda} \
	--recursive \
	--min-qscore ${params.min_qscore} \
	--output-dir . \
  --emit-fastq \
  --emit-summary \
	${params.models} \
	${pod5}
"""
}

process pod5convert{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fast5_folder"
  if (params.basecalling_out != "") {
    publishDir "results/${params.basecalling_out}", mode: 'copy'
  }

  input:
    path(fast5_folder)
  
  output:
    path("converted.pod5"), emit: pod5
  
  script:
"""
pod5 convert fast5 ${fast5_folder}/*.fast5 \
  --output converted.pod5 \
  --threads ${task.cpus} \
  --recursive
"""
}

process demux{
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$fast5_folder"
  if (params.basecalling_out != "") {
    publishDir "results/${params.basecalling_out}", mode: 'copy'
  }

  input:
    path(basecalled)

  output:
    path("*.fastq"), emit: fastq
    path("barcoding_summary.txt"), emit: barcoding_summary

  script:
"""
dorado demux \
  --emit-fastq \
  --emit-summary \
  --threads ${task.cpus} \
  --kit-name ${params.bc_kit} \
  --output-dir . \
  ${basecalled}
"""
}

process basecall {
  container = "${container_url}"
  label "gpus"
  tag "$fast5_folder"
  if (params.basecalling_out != "") {
    publishDir "results/${params.basecalling_out}", mode: 'copy'
  }

  input:
    path(pod5)

  output:
    path("basecalled/"), emit: basecalled

  script:
"""
dorado basecaller \
	-x ${params.cuda} \
	--recursive \
	--min-qscore ${params.min_qscore} \
	--output-dir basecalled/ \
  --kit-name ${params.bc_kit} \
  --no-trim \
	${params.models} \
	${pod5}
"""
}