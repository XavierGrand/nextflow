#!/usr/bin/env Rscript

########################
### 2023-08-21 #########
### H. Polveche ########
########################


################
### Params #####
################

set.seed(123) # Random reproducible

library(optparse)

option_list = list(
  make_option(c("-a", "--cond1"), type="character", default=NULL,
              help="csv file with raw counts for the condition 1", metavar="character"),
  make_option(c("-b", "--cond2"), type="character", default=NULL,
              help="csv file with raw counts for the condition 2 (control condition)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./",
              help="output repertory [default= %default]", metavar="character"),
  make_option(c("-f", "--lfch"), type="double", default=0.4,
              help="log2FoldChange Threshold [default= %default]", metavar="number")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


################
### Packages ###
################

suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(DESeq2, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(viridis, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(gplots, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(openxlsx, quietly = T, warn.conflicts = F))

#################
### Functions ###
################# 

custom.xlsx <- function(wb, tabl, n1, n2, n.sh, namesheet){
  ## custom xlsx final file
  
  addWorksheet(wb, namesheet)
  writeData(wb, sheet = n.sh, tabl, rowNames = FALSE)
  
  hs <- createStyle(fontSize = 10, fgFill = "#dee6ef",
                    halign = "center", valign = "center", textDecoration = "Bold",
                    border = "TopBottomLeftRight") #, textRotation = 45)
  
  co.style <- createStyle(fontSize = 10,
                          halign = "center", valign = "center")
  
  cond1.style <- createStyle(fontSize = 10, fgFill = "#ffffd7",
                             halign = "center", valign = "center")
  cond2.style <- createStyle(fontSize = 10, fgFill = "#e8f2a1",
                             halign = "center", valign = "center")
  sig.style <- createStyle(fontSize = 10, fgFill = "#729fcf",
                           halign = "center", valign = "center")
  
  addStyle(wb, n.sh , co.style, rows = 2:nrow(tabl), 
           cols = 1:ncol(tabl), gridExpand = T)
  addStyle(wb, n.sh , hs, rows = 1, cols = 1:ncol(tabl))
  addStyle(wb, n.sh , hs, rows = 2:nrow(tabl), cols = 1)
  addStyle(wb, n.sh , cond1.style, rows = 2:nrow(tabl), 
           cols = 8:(7+length(n1)), gridExpand = T)
  addStyle(wb, n.sh , cond2.style, rows = 2:nrow(tabl), gridExpand = T,
           cols = (8+length(n1)):(7+length(n1)+length(n2)) )
  addStyle(wb, n.sh, sig.style, rows = 2:nrow(tabl), cols = 6:7, gridExpand = T)
  
  setColWidths(wb, sheet = n.sh, cols = 1, widths = 30)
  setColWidths(wb, sheet = n.sh, cols = 2:ncol(tabl), widths = 10)
  
  return(wb)
}


DEG_conds <- function(cts.cond1 , cts.cond2 , output = "./", LFCh = 0.4 ){
  ### cts.cond1 : rawcount file (csv) to condition (eg. siDDX)
  ### cts.cond2 : rawcounts file (csv) to control consition (eg. siGL2)
  ### output : directory for results
  ### LFCh : Log2FoldChange threshold
  ### DEG with DESeq2 to Find Readthrough
  
  message(" - Formatting input files")
  
  cts.cond1 <- read.csv2(cts.cond1, header = T) %>% 
    select(c(ID, ends_with("_RT")))
  
  names.cond1 <- colnames(cts.cond1)[c(2:ncol(cts.cond1))]

  cts.cond2 <-  read.csv2(cts.cond2, header = T)%>% 
    select(c(ID, ends_with("_RT")))
  
  names.cond2 <- colnames(cts.cond2)[c(2:ncol(cts.cond2))]
  
  cts <- merge(cts.cond1, cts.cond2, by = "ID")
  rownames(cts) <- cts$ID
  cts <- cts[,c(2:ncol(cts))]
  
  coldata <- as.data.frame(matrix(nc=2, nr = length(names.cond1) + length(names.cond2)  ))
  colnames(coldata) <- c("sampleNames", "condition")
  coldata$sampleNames <- c(names.cond1, names.cond2)
  coldata$condition <- as.factor(c(rep("cond1", length(names.cond1)) , rep("cond2", length(names.cond2))    ))
  

  message(" - DEG analysis")  
  cts.1 <- data.frame(cts, moy=NA)
  for (i in 1:nrow(cts.1)){
    cts.1[i, "moy"] <- mean(as.numeric(cts.1[i , c(1:ncol(cts))]))
  }
  cts.2 <- cts.1[which(cts.1$moy > 5),]
  cts.filtre <- cts.2[, c(1:ncol(cts))] 
  colnames(cts.filtre) <- colnames(cts)
  
  cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
  coldata <- coldata[order(coldata$sampleName),]
  
  dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                                design = ~ condition )
  dds <- DESeq(dds)
  
  
  counts_normalise <- counts(dds, norm=T)
  write.table(counts_normalise,
              file=paste0(output, "normcounts.csv"), 
              sep=";", dec=",")
  
  write.table(counts(dds, norm=F),
              file=paste0(output, "rawcounts.csv"), 
              sep=";", dec=",")
  
  resGA <- results(dds, contrast=c("condition","cond1", "cond2"), 
                   lfcThreshold=LFCh, altHypothesis="greaterAbs")  # lfcThreshold=0.4
  message("   \n")
  message (paste0("|Log2FoldChange| Threshold > ", LFCh))
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
  
  write.table(resGA, paste0(output, "cond1-vs-cond2_all_genes.csv"), 
              sep =  ";", dec = ",")
  write.table(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),], 
              paste0(output, "cond1-vs-cond2_sig_log2FC", LFCh,"_padj005_BM20.csv"),
              sep = ";", dec = ",")
  
  
  message(" - Final Excel File")
  
  resGA.m <- merge(as.data.frame(resGA), round(counts_normalise, digit = 2), by = "row.names")
  colnames(resGA.m)[1] <- "ID"
  resGA.m[,c(2:5)] <- round(resGA.m[,c(2:5)], digits = 2)
  res.sig.m <- resGA.m[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]
    
  wb <- createWorkbook()
  wb <- custom.xlsx(wb=wb, tabl = res.sig.m, n1 = names.cond1, 
                    n2 = names.cond2, n.sh = 1, namesheet = "sig")
  wb <- custom.xlsx(wb, resGA.m, names.cond1, names.cond2, 2, "all")
  saveWorkbook(wb, paste0(output,"cond1-vs-cond2_DESeq2_results.xlsx"), overwrite = TRUE)
  
  message(" - Done")
}



#################
### __main__ ####
################# 


# cts.cond1 <- "./results/siDDX/rawcounts.csv"
# cts.cond2 <- "./results/siGL2/rawcounts.csv"
# output <- "./results/siDDX-vs-siGL2/"
# LFCh <- 0.4


DEG_conds(cts.cond1 = opt$cond1, cts.cond2 = opt$cond2, output = opt$out, LFCh = opt$lfch)
# DEG_conds(cts.cond1 , cts.cond2 , output, LFCh)









