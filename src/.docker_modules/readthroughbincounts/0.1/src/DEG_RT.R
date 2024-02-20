#!/usr/bin/env Rscript

########################
### 2023-08-03 #########
### H. Polveche ########
########################


################
### Params #####
################

set.seed(123) # Random reproducible



library(optparse)


option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Directory with kallisto folders", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./",
              help="output repertory [default= %default]", metavar="character"),
  make_option(c("-f", "--lfch"), type="double", default=4,
              help="log2FoldChange Threshold [default= %default]", metavar="number"),
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="Number of process parallelization [default= %default]")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input directory).n", call.=FALSE)
}

n.worked <- opt$threads



################
### Packages ###
################

suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(doFuture, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(DESeq2, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(viridis, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(gplots, quietly = T, warn.conflicts = F))



#################
### Functions ###
################# 

registerDoFuture()
plan(cluster, workers = n.worked)

format_table_RT_1file <- function(file = "./abundance.tsv", name = "toto" ){
  ### file abundance kallisto output
  ### passer de 1 colonne a 2 colonne gene et gene_RT
  
  kal <- read.csv2(file, header = T, sep = "\t", dec = ".") %>% 
    select(target_id, est_counts) %>% 
    mutate( RTdetect = str_detect(string = target_id, pattern = "_RT$"))
  
  kal.genes <- kal %>% 
    filter( RTdetect == FALSE) %>% 
    select(target_id, est_counts)
  colnames(kal.genes)[2] <- name
  
  kal.RT <- kal %>% 
    filter( RTdetect == TRUE) %>% 
    select(target_id, est_counts) %>% 
    mutate(target_id = str_replace(target_id, "_RT", "") )
  colnames(kal.RT)[2] <- paste0(name,"_RT")
  
  res <- merge(kal.genes, kal.RT, by = "target_id")
  
  return(res)
}



format_tables_toDESeq2 <- function(dir = "./"){
  ### dir : directory with kallisto folders
  
  if ( !(dir.exists(dir)) ){
    stop(paste0("The directory '", dir,"' does not exist."))
  }
  
  list_dir <- list.files(path = dir, full.names = F)
  
  dfs <- foreach(i = 1:length(list_dir), .combine=cbind)  %dopar% #, .combine=cbind)
    format_table_RT_1file(file = paste0(dir, list_dir[i], "/abundance.tsv"), name = list_dir[i] )
  
  colnames(dfs)[1] <- "ID"
  dfs <- dfs %>% 
    select(!target_id)
  
  dfs[ ,c(2:ncol(dfs))] <- round( dfs[ ,c(2:ncol(dfs))] )
  
  #reg.ex <- commun_expression(list_dir)
  coldata <- data.frame( sampleName = colnames(dfs)[c(2:ncol(dfs))], condition = rep(c("canonical", "RT"), (ncol(dfs) - 1)/2))
  rownames(coldata) <- coldata$sampleName
  
  return(list("counts" = dfs, "coldata" = coldata)) 
}



DEG_RT <- function(input = "./data/kallisto/", output = "./results/", LFCh = 4 ){
  ### input : directory with kallisto folders
  ### output : directory for results
  ### LFCh : Log2FoldChange threshold
  ### DEG with DESeq2 to Find Readthrough
  
  message(" - Formatting input files")
  mat <- format_tables_toDESeq2(dir = input)
  
  cts <- mat$counts
  write_excel_csv(cts, paste0(output, "rawcounts.csv"), col_names = T, delim = ";")
  
  rownames(cts) <- cts$ID
  cts <- cts[,c(2:ncol(cts))]
  
  coldata <- mat$coldata
  write_excel_csv(coldata, paste0(output, "coldata.csv"), col_names = T, delim = ";")
  coldata$condition <- as.factor(coldata$condition)

  message(" - DEG analysis")  
  cts.1 <- data.frame(cts, moy=NA)
  for (i in 1:nrow(cts.1)){
    cts.1[i, "moy"] <- mean(as.numeric(cts.1[i , c(1:ncol(cts))]))
  }
  cts.2 <- cts.1[which(cts.1$moy > 5),]
  cts.filtre <- cts.2[, c(1:ncol(cts))] # 12949
  colnames(cts.filtre) <- colnames(cts)
  
  colnames(cts.filtre) <- c(paste0("X", colnames(cts.filtre)))
  rownames(coldata) <- c(paste0("X", rownames(coldata)))
  
  cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
  coldata <- coldata[order(coldata$sampleName),]
  
  dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                                design = ~ condition )
  dds <- DESeq(dds)

  
  counts_normalise <- counts(dds, norm=T)
  write.table(counts_normalise,
              file=paste0(output, "normcounts.csv"), 
              sep=";", dec=",")
  
  

  resGA <- results(dds, contrast=c("condition","RT", "canonical"), 
                   lfcThreshold=LFCh, altHypothesis="greater")  # lfcThreshold=0.4
  message("   \n")
  message (paste0("|Log2FoldChange| Threshold > 4"))
  message(paste0( "Number of genes with Readthrough : ",
    nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]) # Up
  ))
  
  
  message(" - Plots")
  rld <- rlogTransformation(dds, blind=TRUE, fitType = "local")
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE, fitType = "local")
  vstMat = assay(vsd)
  
  pdf(paste0(output, "/DEG_plots.pdf"))
  
    condcols=c("firebrick","steelblue")
    names(condcols)=unique(coldata$condition)
    
    ylim <- c(-10,10)
    drawLines <- function() abline(h=c(-0.4,0.4),col="steelblue",lwd=2)
    plotMA(resGA, ylim=ylim); drawLines()
    
    hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])
    
    plot(resGA$log2FoldChange,-log(resGA$padj,10),
         ylab='-log10(Adjusted P)',
         xlab="Log2 FoldChange",
         pch=19,cex=0.5, col = "dimgray"
    )      
    
    points(resGA[hits,'log2FoldChange'],
           -log(resGA[hits,'padj'],10),
           pch=19,
           cex=0.5,
           col="steelblue"
    )
    abline(h=-log10(0.05),lty=3)
    abline(v=-LFCh,lty=3)
    abline(v=LFCh,lty=3)
    
    
    barplot(colSums(counts(dds, normalized=F)), col=condcols[as.factor(coldata$condition)], 
            las=2,cex.names=0.4,
            main='Pre Normalised Counts')
    
    barplot(colSums(counts(dds, normalized=T)), col=condcols[as.factor(coldata$condition)], 
            las=2,cex.names=0.4,
            main='Post Normalised Counts')
    
    pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    ggplot(pcaData, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
    
    sampleDists <- dist( t( assay(rld) ) )
    sampleDistMatrix <- as.matrix( sampleDists )
    colours <- inferno(255)
    heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15), 
                     cexRow=0.5, cexCol=0.5)
    
    plotDispEsts(dds)
  
  dev.off()
  
  write.table(resGA, paste0(output, "RT-vs-Canonical_all_genes.csv"), 
                  sep =  ";", dec = ",")
  write.table(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),], 
                  paste0(output, "RT-vs-Canonical_sig_log2FC", LFCh,"_padj005_BM20.csv"),
                  sep = ";", dec = ",")
    
  message(" - Done")
}



#################
### __main__ ####
################# 

DEG_RT(input = opt$dir, output = opt$out, LFCh = opt$lfch)
#DEG_RT(input = "./data/kallisto/siGL2/", output = "./results/", LFCh = 4)




