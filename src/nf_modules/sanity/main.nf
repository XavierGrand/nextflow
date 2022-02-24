container_url="mlepetit/sanity:latest"

params.sanity_out=""
params.sanity=""

process normalization_sanity
        {

        container="${container_url}"
        label  "big_mem_multi_cpus"
        if (params.sanity_out != "") {
		publishDir "results/${params.sanity_out}", mode: 'copy'

	}
else {
          publishDir "results/normalize_matrix/", mode: 'copy'

           }

	input:

               tuple val(id_mtx), path(raw_filtered_mtx)   
               

        output:

               tuple val(id_mtx),path("log_transcription_quotients.txt"), emit: normalize_filtered_mtx
               tuple val(id_mtx), path("ltq_error_bars.txt")  ,emit: ltq_error

        script:

        """
        Sanity -f ${raw_filtered_mtx} -n ${task.cpus} ${params.sanity} 
        """
        }
