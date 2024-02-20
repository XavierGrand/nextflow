#!/usr/bin/env Rscript

########################
### 2023-07-27 #########
### H. Polveche ########
########################

### script --vanilla ./src/readthrough_bin.R


#library("optparse")

# option_list = list(
#   make_option(c("-f", "--file"), type="character", default=NULL, 
#               help="dataset file name", metavar="character"),
#   make_option(c("-o", "--out"), type="character", default="out.txt", 
#               help="output file name [default= %default]", metavar="character")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 

# if (is.null(opt$file)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }

# ## program...
# df = read.table(opt$file, header=TRUE)
# num_vars = which(sapply(df, class)=="numeric")
# df_out = df[ ,num_vars]
# write.table(df_out, file=opt$out, row.names=FALSE)

#Rscript --vanilla yasrs.R -f iris.txt -o out.txt


#################
### Variables ###
#################

fileIN <- "./data/readthrough_range.bed"
#pathOUT <- "./results/"
#n.worked <- 6
#suffix <- "toto"



################
### Packages ###
################

library(tidyverse)
library(progress) # progress bar
library(doFuture) # parallelism
library(GenomicRanges)
library(valr)
library(Rsamtools)
library(bamsignals)
library(rtracklayer)


##############
### Config ###
##############


set.seed(123) # Random reproducible

## registerDoFuture()
## plan(cluster, workers = n.worked)
#plan(multisession, workers = n.worked)

source(file = "./src/readthrough_bin.R")

# message("Number of parallel workers: ", nbrOfWorkers())

bed <- read.csv2(fileIN, sep ="\t", header = F)
### LAAAA
vars <- normalize.stranded(bedVar = bed, suf = "tutu", pathOUT = "./results/")
df.bin <- bintab(vars = vars, bin = 50, n.worked = 6, suf = "tutu", pathOUT = "./results/")

cp <- countBAMwtBED(bedfile = "./results/Regions_50bins_tutu.bed",
                        bamPath = "./data/BAM/5Y_siDDX5-17_B1/5Y_siDDX5-17_B1.bam")
counts <- cp$counts
coverage <- cp$profile
counts.3 <- counts[which(counts$moy >= 3),]







####################
### Progress Bar ###
####################

n_iter <- 200
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = n_iter,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 200)      # Width of the progress bar

for(i in 1:n_iter) {
  
  # Updates the current state
  pb$tick()
  
  #---------------------
  # Code to be executed
  #---------------------
  
  Sys.sleep(0.1) # Remove this line and add your code
  
  #---------------------
  
}



