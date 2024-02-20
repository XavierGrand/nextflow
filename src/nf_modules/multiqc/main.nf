// multiqc generate nice html report combining lots of differents bioinformatics
// tools report.
// 
// EXAMPLE:

/*
include { multiqc } 
  from './nf_modules/multiqc/main'
  addParams(
    multiqc_out: "QC/"
  )

multiqc(
  report_a
  .mix(
    report_b,
    report_c,
    report_d
  )
)
*/

version = "1.11"
container_url = "lbmc/multiqc:${version}"

params.multiqc = ""
params.multiqc_out = "QC/"
workflow multiqc {
  take:
    report
  main:
    report
    .map{it ->
      if (it instanceof List){
        if(it.size() > 1) {
          it[1]
        } else {
          it[0]
        }
      } else {
        it
      }
    }
    .unique()
    .flatten()
    .set { report_cleaned }
    multiqc_default(report_cleaned.collect())

  emit:
  report = multiqc_default.out.report
}

process multiqc_default {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  if (params.multiqc_out != "") {
    publishDir "results/${params.multiqc_out}", mode: 'copy'
  }

  input:
    path report 

  output:
    path "*multiqc_*", emit: report

  script:
"""
multiqc ${params.multiqc} -f .
"""
}
