version = "5.0.11"
container_url = "lbmc/guppy-gpu:${version}"

params.basecalling_out = ""
params.flowcell = ""
params.kit = ""
params.gpu_runners_per_device = 16
process basecall_fast5 {
  container = "${container_url}"
  // Need to create a profile using GPUs
  label ""
  tag "$file_id"
  if (params.basecalling_out != "") {
    publishDir "results/${params.basecalling_out}", mode: 'copy'
  }

  if (params.flowcell == "") {
      errorFlowcell << "WARNING ! No Flowcell type given..."
      errorFlowcell.view()
  }

  if (params.kit == "") {
      errorKit "WARNING ! No kit type given..."
      errorKit.view()
  }

  input:
    tuple val(file_id), path(fast5)

  output:
    tuple val(file_id), path("*.fastq*"), emit: fastq

  script:
"""
guppy_basecaller --compress_fastq -x "cuda:all" --min_qscore 7.0 \
    -i ${path(fast5)} \
    -s ${params.basecalling_out} \
    --gpu_runners_per_device ${params.gpu_runners_per_device} \
    --flowcell ${params.flowcell} \
    --kit ${params.kit}
"""
}