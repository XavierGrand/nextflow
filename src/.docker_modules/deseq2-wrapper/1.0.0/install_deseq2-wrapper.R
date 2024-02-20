#!/usr/bin/env Rscript
install.packages("devtools")
library(devtools)
install_gitlab("LBMC/regards/deseq2-wrapper", host = "https://gitbio.ens-lyon.fr", quiet = FALSE)
