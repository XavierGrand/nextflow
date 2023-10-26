#!/usr/bin/env Rscript

########################
### 2023-08-03 #########
### H. Polveche ########
########################


################
### Params #####
################

library(optparse)

option_list = list(
   make_option(c("-f", "--file"), type="character", default=NULL,
               help="GTF File", metavar="character"),
   make_option(c("-o", "--out"), type="character", default="./",
               help="output repertory [default= %default]", metavar="character"),
   make_option(c("-s", "--suffix"), type="character", default="readtrhough",
               help="A suffix to specify your results [default= %default]", metavar="character"),
   make_option(c("-r", "--RT"), type="integer", default=10000,
               help="Size of Readthrough search zone after last gene exon (#bases) [default= %default]", metavar="number"),
   make_option(c("-c", "--chr"), type="character", default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,X,Y,MT",
               help="List of chromosomes you want to keep (default: humain chr) [default= %default]")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$file)){
   print_help(opt_parser)
   stop("At least one argument must be supplied (input GTF file).n", call.=FALSE)
}


################
### Packages ###
################

suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(GenomicRanges, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(valr, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(rtracklayer, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(progress, quietly = T, warn.conflicts = F) )

#################
### Functions ###
################# 

gtfConvert <- function(file, out , suf, RT, chr){
  ### file gtf with grep function
  ### e.g : 
  ### grep -P "\tgene\t" Homo_sapiens.GRCh37.75.gtf > Homo_sapiens.GRCh37.75.genes.gtf
  
  vec.chr <- str_split(chr, ",")
  vec.chr <- vec.chr[[1]]
  
  message("Loading GTF file ( 1 / 6 )")
  gtf <- rtracklayer::import(file)
  gtf <- gtf[which(gtf$gene_biotype %in% "protein_coding"),]
  gtf$score <- 0
  gtf$name <- paste0(gtf$gene_name,"_", gtf$gene_id)
  
  
  gtf10kb <- flank(gtf, start = FALSE, width = RT)
  gtf10kb$name <- paste0(gtf10kb$name, "_RT")
  
  message("Export classical Bed file ( 2 / 6 )")
  
  rtracklayer::export(gtf, paste0(out, "Genes_", RT ,"bases_", suf ,".bed"), format="bed") #,object = gtf, genome = "hg19", format = "bed15", name = "GRCh37.75_genes")
  rtracklayer::export(gtf10kb, paste0(out, "GenesAndRT_", RT ,"bases_", suf ,".bed"), format="bed")
  
  gtf.r <- read_bed(paste0(out, "Genes_", RT ,"bases_", suf ,".bed")) %>% 
    filter(chrom %in% vec.chr)
  
  # On ne conserve que les chr canoniques et on place à 0 si strand (-) au lieu de valeur neg
  gtf10kb.r <- read_bed(paste0(out, "GenesAndRT_", RT ,"bases_", suf ,".bed")) %>% 
    filter(chrom %in% vec.chr) %>% 
    mutate(start = if_else(start < 0 , 0, start)) 
  file.remove(paste0(out, "GenesAndRT_", RT ,"bases_", suf ,".bed")) 
  
  gtf10kb.r.val <- as.data.frame(gtf10kb.r) %>% 
    rename(X4 = "symbol" ,
           X6 = "strand")
  gtf.r.val <- as.data.frame(gtf.r) %>% 
    rename(X4 = "symbol" ,
           X6 = "strand")
  
  message(paste0("Create ", RT,"kb extension ( 3 / 6 )"))
  
  # intersection pour avoir les zones RT chevauchant avec des genes
  overl <- valr::bed_intersect(gtf10kb.r.val, gtf.r.val, suffix = c("_10kb", "")) %>% 
    filter(.overlap > 0) %>% 
    filter((strand_10kb == strand)) %>% 
    mutate(end_10kb = if_else(strand_10kb == "+", start, end_10kb)) %>% 
    mutate(start_10kb = if_else(strand_10kb == "-", end, start_10kb)) %>% 
    select(chrom, start_10kb, end_10kb, symbol_10kb, X5_10kb, strand_10kb) %>% 
    filter(start_10kb - end_10kb < 0)
  
  message("Overlap on genes : coordinate correction ( 4 / 6 )")
  
  # on corrige pour que la fin du RT ne soit pas de 10kb mais avant le debut du gene n+1
  gtf10kb.r.val2 <- gtf10kb.r.val 
  pb <- progress_bar$new(format = " [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]", #(:spin)
                         total = nrow(overl),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  for (i in 1:nrow(overl)){
    pb$tick()
    if (overl[i, "symbol_10kb"] %in% gtf10kb.r.val2[,"symbol"]){
      gtf10kb.r.val2[which(gtf10kb.r.val2$symbol %in% 
                             overl[i, "symbol_10kb"]), c("start", "end")] <- overl[i, c("start_10kb", "end_10kb")] 
    }
  }
  
  message("Overlap on genes : Stranded correction ( 5 / 6 )")
  pb <- progress_bar$new(format = " [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]", # (:spin) 
                         total = nrow(gtf10kb.r.val2),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  # On garde les coordonnees des RT uniquement
  write_excel_csv(gtf10kb.r.val2, paste0(out, "RT-Only_", RT ,"bases_", suf ,".bed"), 
                  delim = "\t", col_names = F, quote="none")
  
  # On met le debut du RT au meme niveau que le gene associe
  for (i in 1:nrow(gtf10kb.r.val2)){
    pb$tick()
    if(gtf10kb.r.val2[i, "strand"] == "+"){
      gtf10kb.r.val2[i, "start"] <- gtf.r.val[which(gtf.r.val$symbol %in% str_replace(gtf10kb.r.val2[i, "symbol"], "_RT", "")), "start"] #which(gtf.r.val$symbol %in% str_replace(gtf10kb.r.val2[i, "symbol"], "_RT", ""))
      
    } else {
      gtf10kb.r.val2[i, "end"] <- gtf.r.val[which(gtf.r.val$symbol %in% str_replace(gtf10kb.r.val2[i, "symbol"], "_RT", "")), "end"] #which(gtf.r.val$symbol %in% str_replace(gtf10kb.r.val2[i, "symbol"], "_RT", ""))
      
    }
    
    
  }
  
  message("Export Extension Bed file ( 6 / 6 )")

  gtf10kb.r.val2 <- rbind(gtf10kb.r.val2, gtf.r.val)
  write_excel_csv(gtf10kb.r.val2, paste0(out, "GenesAndRT_", RT ,"bases_", suf ,".bed"), 
                  delim = "\t", col_names = F, quote="none")  
  
  #return(list("Genes" = gtf.r.val, "RT" = gtf10kb.r.val2))
}


#################
### __main__ ####
################# 

gtfC <- gtfConvert(file =  opt$file ,
                   out  =  opt$out ,
                   suf  =  opt$suffix ,
                   RT   =  opt$RT , 
                   chr  =  opt$chr
)


