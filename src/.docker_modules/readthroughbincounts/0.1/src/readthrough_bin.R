########################
### 2023-07-24 #########
### H. Polveche ########
########################


#################
### Functions ###
#################



normalize.stranded <- function(bedVar, suf = "suf", pathOUT = "./"){
  ### bedVar : tableau de sortie du programme python ReadThrough X. Grand
  ### https://gitbio.ens-lyon.fr/xgrand/readthrough/-/tree/master/
  ### suf : suffix to output bed
  ### pathOUT : chemin pour les resultats 
  ### Formatage du tableau de sortie du programme de X. Grand
  
  
  vec.chr <- c(1:23, "X", "Y", "MT")
  
  dt <- data.frame(bedVar, start.exon=NA, end.exon=NA, region500=NA, 
                   end.rt = NA, V14=NA, V15=NA)
  dt <- dt %>% 
    mutate(start.exon = case_when(V6 == "+" ~ V9,
                                  V6 == "-" ~ V10
                                  )) %>% 
    mutate(end.exon = case_when(V6 == "+" ~ V10,
                                V6 == "-" ~ V9
                                )) %>% 
    mutate(region500 = case_when(V6 == "+" ~ V2,
                                V6 == "-" ~ V3
                                )) %>%
    mutate(end.rt = case_when(V6 == "+" ~ V3,
                              V6 == "-" ~ V2
                                )) %>% 
    mutate(V14 = case_when(V6 == "+" ~ V2,
                          V6 == "-" ~ V3
                                )) %>%
    mutate(V15 = case_when(V6 == "+" ~ V3,
                          V6 == "-" ~ V2
                                )) %>% 
    select(V1, V4, V6, V7, start.exon, end.exon, region500, end.rt, V14, V15) %>% 
    dplyr::rename("chr" = "V1",
           "strand" = "V6",
           "gene_symbol" = "V7",
           "ENSG" = "V4") %>% 
    filter(chr %in% vec.chr)
  
  bedT <- dt[,c("chr", "V14", "V15", "gene_symbol")]
  colnames(bedT) <- c("chr","start","end",  "gene_symbol")
  dt <- dt %>% 
    select(chr, strand, gene_symbol, ENSG, start.exon, end.exon, region500, end.rt)
  write_tsv(bedT, paste0(pathOUT,"Regions10500bases_toBedtools_", suf, ".bed"), 
            col_names = F )
  write_excel_csv(dt, paste0(pathOUT,"Regions10500bases_Coordinates_infos_", suf, ".csv"), 
                  delim = ";" )
 
  return(list("table" = dt, "bed" = bedT))
}



OneGeneBin <- function(bedT, tab, bin = 50){
  ### bedT : sortie 'bed' de la fonction normalize.stranded(), 
  ###     ligne i de la fonction bintab()
  ### tab : sortie 'table' de la fonction normalize.stranded(), 
  ###     ligne i de la fonction bintab()
  ### bon : taille d'un bin pour le decoupage
  ### Pour un ligne ( = chaque gene ) on va couper la région 500+10k bases
  ###     en bin de {bin}
  
  if (  ((abs(bedT[1, "start"] - bedT[1, "end"])) %% bin) %in% 0){
    stra <- tab[1, "strand"] 
    end.ex <- tab[1, "end.exon"]
    repet <- abs(bedT[1, "start"] - bedT[1, "end"]) / bin + 1
    tabi <- as.data.frame(matrix(nc=11, nr = repet -1))
    colnames(tabi) <- c("chr", "start", "end", "ID", "bins", "strand", 
                        "gene_symbol", "ENSG", "category","s.bed", "e.bed")
    tabi[, "gene_symbol"] <- rep(tab[1, "gene_symbol"], repet -1)
    tabi[, "chr"] <- rep(tab[1, "chr"], repet -1)
    tabi[, "ENSG"] <- rep(tab[1, "ENSG"], repet -1)
    tabi[, "strand"] <- rep(stra, repet -1)
    
    
    if(stra == "+"){
      bin2 <- bin
      
      vec <- seq(bedT[1, "start"], bedT[1, "end"], by = bin2)
      tabi[, "start"] <- vec[1:(repet-1)]
      tabi[, "s.bed"] <- vec[1:(repet-1)]
      tabi[, "end"] <- vec[2:repet] - 1
      tabi[, "e.bed"] <- vec[2:repet] - 1
      tabi[, "bins"] <- seq(1:(repet-1))
      
    } else {
      bin2 <- (-1) * bin
      vec <- seq(bedT[1, "start"], bedT[1, "end"], by = bin2)
      tabi[, "start"] <- vec[1:(repet-1)]
      tabi[, "end"] <- vec[2:repet] - 1
      tabi[, "e.bed"] <- vec[1:(repet-1)]
      tabi[, "s.bed"] <- vec[2:repet] - 1
      tabi[, "bins"] <- seq(1:(repet-1))
      
    }
   
    
    tabi <- tabi %>% 
      mutate(ID = paste0(gene_symbol, "_", bins)) %>% 
      mutate(category = case_when(stra == "+" ~ ifelse(end.ex >= start, "EXON", "RT"),
                                  stra == "-" ~ ifelse(end.ex >= start, "RT", "EXON")))
    
  }
}




bintab <- function(vars, bin = 50, n.worked = 6, suf = "suf", pathOUT = "./"){
  ### binage de tous les genes , parralelisation et optimisation de temps de calcul
  ### vars : result fonction normalize.stranded()
  ### suf : suffix to output bed
  ### pathOUT : chemin pour les resultats 
  
  tab <- vars$table
  bedT <- vars$bed
  
  set.seed(123) # Random reproducible
  
  registerDoFuture()
  plan(cluster, workers = n.worked)

  tabo <- foreach(i = 1:nrow(tab)) %dopar% { 
    OneGeneBin(bedT[i, ], tab[i, ], bin = 50)
  } 

  tabo <- compact(tabo)
  df <- dplyr::bind_rows(tabo)
  
  write_tsv(df[,c("chr", "s.bed", "e.bed", "ID")], paste0(pathOUT,"Regions_", bin, "bins_", 
                         suf, ".bed"), col_names = F )
  write_excel_csv(df, paste0(pathOUT,"Regions_", bin, "bins_", suf, ".csv"), delim = ";" )

 return(df)

}


countBAMwtBED <- function(bedfile, bamPath){
  ### bedfile : result de bintab()
  ### "./results/Regions_50bins_tutu.bed"
  ### bamPath : chemin du fichier bam
  ### "./results/DDX_6_aligned_sorted.filter.bin50.bam"
  
  gr_obj <-  rtracklayer::import(bedfile)
  bamFile <- Rsamtools::BamFile(bamPath)
  
  sigs.c <- bamCount(bamPath, gr_obj, mapq = 10, verbose=FALSE)
  sigs.p <- bamCoverage(bamPath, gr_obj, mapq = 10, verbose=FALSE)
  
  df.bin.c <- data.frame(df.bin, count=NA, moy=NA)
  df.bin.c[, "count"] <- sigs.c 
  df.bin.c[, "moy"] <- round(sigs.c / 50 )
  

  return(list("counts" = df.bin.c, "coverage" = sigs.p))
}

