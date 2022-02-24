version = "3.0.0"
container_url = "mlepetit/cellphonedb:latest"

params.cellphonedb = ""
params.cellphonedb_out = ""
params.pval=""
params.thres=""
params.iterations=""
params.gene_id=""



workflow cellphone_statistical_analysis {
  take:
    meta
    counts

  main:

cpdb_methods_stats(meta,counts)
cpdb_plot_dot_plot(cpdb_methods_stats.out.means,cpdb_methods_stats.out.pvalues)
cpdb_plot_heatmap(cpdb_methods_stats.out.pvalues)


  emit:
    means = cpdb_methods_stats.out.means
    pvalues = cpdb_methods_stats.out.pvalues
    deconvoluted = cpdb_methods_stats.out.deconvoluted
    significant_means = cpdb_methods_stats.out.significant_means
    dot_plot = cpdb_plot_dot_plot.out.dot_plot
    heatmap = cpdb_plot_heatmap.out.heatmap
    heatmap_log = cpdb_plot_heatmap.out.heatmap_log
    count_network = cpdb_plot_heatmap.out.count_network
    interactions_count = cpdb_plot_heatmap.out.interactions_count


}










process cpdb_methods_stats {
  container = "${container_url}"
  label "big_mem_multi_cpus"
    if (params.cellphonedb_out != "") {
    publishDir "results/${params.cellphonedb_out}", mode: 'copy'
  }

  input:
    tuple val(id_mtx), path(meta)
    tuple val(id_mtx), path(counts)
  
  output:
   tuple val(id_mtx), path("out/means.txt"), emit: means
   tuple val(id_mtx), path("out/pvalues.txt"), emit: pvalues
   tuple val(id_mtx), path("out/deconvoluted.txt"), emit: deconvoluted
   tuple val(id_mtx), path("out/significant_means.txt"), emit: significant_means
  
script:
  """
cellphonedb method statistical_analysis ${params.meta} ${params.counts} --counts-data ${params.gene_id}  --threads ${task.cpus} --iterations ${params.iterations} --pvalue ${params.pval} --threshold ${params.thres}

  """
}


process cpdb_plot_dot_plot {
  container = "${container_url}"
  label "big_mem_mono_cpus"
    if (params.cellphonedb_out != "") {
    publishDir "results/${params.cellphonedb_out}", mode: 'copy'
  }

  input:
    tuple val(id_mtx), path(means)
    tuple val(id_mtx), path(pvalues)

  output:
   tuple val(id_mtx), path("out/plot.pdf"), emit: dot_plot

script:
  """
mkdir ./out
cellphonedb plot dot_plot --means-path ${means} --pvalues-path ${pvalues} 

  """
}
 
process cpdb_plot_heatmap {
  container = "${container_url}"
  label "big_mem_multi_cpus"
    if (params.cellphonedb_out != "") {
    publishDir "results/${params.cellphonedb_out}", mode: 'copy'
  }

  input:
    tuple val(id_mtx), path(pvalues)

  output:
   tuple val(id_mtx), path("out/heatmap_count.pdf"), emit: heatmap
 tuple val(id_mtx), path("out/heatmap_log_count.pdf"), emit: heatmap_log
 tuple val(id_mtx), path("out/count_network.txt"), emit: count_network
 tuple val(id_mtx), path("out/interaction_count.txt"), emit: interactions_count

script:
 
 """
mkdir ./out
cellphonedb plot heatmap_plot --pvalues-path ${pvalues} --pvalue ${params.pval} ${params.meta}

  """
}
