#!/usr/bin/env Rscript

########################
### 2023-08-04 #########
### H. Polveche ########
########################


################
### Params #####
################

set.seed(123) # Random reproducible


# library(optparse)
# 
# 
# option_list = list(
#   make_option(c("-d", "--dir"), type="character", default=NULL,
#               help="Directory with kallisto folders", metavar="character"),
#   make_option(c("-o", "--out"), type="character", default="./",
#               help="output repertory [default= %default]", metavar="character"),

# );
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# 
# if (is.null(opt$dir)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input directory).n", call.=FALSE)
# }


################
### Packages ###
################

suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(rtracklayer, quietly = T, warn.conflicts = F))
#suppressPackageStartupMessages(library(doFuture, quietly = T, warn.conflicts = F))

library(bamsignals)


rangeBed <- function(resGA = "./results/siDDX5-17/RT-vs-Canonical_sig_log2FC4_padj005_BM20.csv", 
                     bed = "./results/RT-Only_10000bases_essai.bed"){
  ### ne garder que les coordonnees de RT qui sont DEG sig 
  
  
  
  
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


