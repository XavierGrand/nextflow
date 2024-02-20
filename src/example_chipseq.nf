
nextflow.enable.dsl=2

include {
  fastp
} from './nf_modules/fastp/main'

workflow csv_parsing {
  if (params.csv_path.size() > 0) {
    log.info "loading local csv files"
    Channel
      .fromPath(params.csv_path, checkIfExists: true)
      .ifEmpty { error 
      log.error """
    =============================================================
      WARNING! No csv input file precised.
      Use '--csv_path <file.csv>'
      Or '--help' for more informations
    =============================================================
    """
      }
      .splitCsv(header: true, sep: ";", strip: true)
      .flatMap{
        it -> [
          [(it.IP + it.WCE).md5(), "IP", "w", file(it.IP)],
          [(it.IP + it.WCE).md5(), "WCE", "w", file(it.WCE)]
        ]
      }
      .map{ it ->
        if (it[1] instanceof List){
          it
        } else {
          [it[0], [it[1]], it[2], it[3], [it[4]]]
        }
      }
      .map{
        it ->
        if (it[1].size() == 2){ // if data are paired_end
          [
            "index": it[0],
            "group": ref_order(it),
            "ip": it[2],
            "type": it[3],
            "id": read_order(it)[0],
            "file": read_order(it)
          ]
        } else {
          [
            "index": it[0],
            "group": it[1].simpleName,
            "ip": it[2],
            "type": it[3],
            "id": it[4].simpleName,
            "file": [it[4].simpleName, it[4]]
          ]
        }
      }
      .set{input_csv}
  } else {
    log.info "loading remotes SRA csv files"
    Channel
      .fromPath(params.csv_sra, checkIfExists: true)
      .ifEmpty { error 
      log.error """
    =============================================================
      WARNING! No csv input file precised.
      Use '--csv_path <file.csv>' or
      Use '--csv_SRA <file.csv>'
      Or '--help' for more informations
    =============================================================
    """
      }
      .splitCsv(header: true, sep: ";", strip: true)
      .flatMap{
        it -> [
          [[it.IP_w + it.WCE_w + it.IP_m + it.WCE_m], t.IP_w, "IP", "w", it.IP_w],
          [[it.IP_w + it.WCE_w + it.IP_m + it.WCE_m], it.IP_w, "WCE", "w", it.WCE_w],
          [[it.IP_w + it.WCE_w + it.IP_m + it.WCE_m], it.IP_w, "IP", "m", it.IP_m],
          [[it.IP_w + it.WCE_w + it.IP_m + it.WCE_m], it.IP_w, "WCE", "m", it.WCE_m]
        ]
      }
      .map{
        it ->
        if (it[1].size() == 2){ // if data are paired_end
          [
            "index": (
              it[0][0][0].simpleName +
              it[0][0][1].simpleName +
              it[0][0][2].simpleName +
              it[0][0][3].simpleName
            ).md5(),
            "group": it[1][0].simpleName,
            "ip": it[2],
            "type": it[3],
            "id": it[4][0].simpleName[0..-4],
            "file": [it[4][0].simpleName[0..-4], it[4]]
          ]
        } else {
          [
            "index": (
              it[0][0].simpleName +
              it[0][1].simpleName +
              it[0][2].simpleName +
              it[0][3].simpleName
            ).md5(),
            "group": it[1].simpleName,
            "ip": it[2],
            "type": it[3],
            "id": it[4].simpleName,
            "file": [it[4].simpleName, it[4]]
          ]
        }
      }
      .set{input_csv}
  }
  emit:
  input_csv
}


workflow {

}