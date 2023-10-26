#!/bin/Rscript

### Author: Xavier Grand
### Date: 23/10/2023

# Packages installation:

list.of.packages <- c("BiocManager", "progress", "doFuture", "viridis",
                      "gplots", "openxlsx", "optparse")

new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T)

list.of.BiocPackages <- c("rtracklayer", "GenomicRanges", "DESeq2", "Rsamtools", "bamsignals")
new.BiocPackages <- list.of.BiocPackages[!(list.of.BiocPackages %in% installed.packages()[,"Package"])]

if (length(new.BiocPackages) > 0) {
  BiocManager::install(new.BiocPackages, dependencies = TRUE, ask = FALSE)
}

install.packages("valr", dependencies = T)
