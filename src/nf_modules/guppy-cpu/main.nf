version = "5.0.11"
container_url = "lbmc/guppy-cpu:${version}"

params.basecalling_out = ""
params.flowcell = "FLO-MIN106"
params.kit = "SQK-PCS109"
params.cpu_threads_per_caller = 4
params.num_callers = 1
process basecall_fast5 {
  container = "${container_url}"
  label "big_mem_multi_cpus"
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
guppy_basecaller --compress_fastq \
    -i ${path(fast5)} \
    -s ${params.basecalling_out} \
    --cpu_threads_per_caller ${params.cpu_threads_per_caller} \
    --num_callers ${params.num_callers} \
    --flowcell ${params.flowcell} \
    --kit ${params.kit}
"""
}