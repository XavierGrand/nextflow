#!/usr/bin/env Rscript
list.of.packages <- c("optparse", "dplyr", "tidyr", "stringr", "BiocManager",
                      "tidyverse", "tibble", "purrr", "furrr", "future", "vsn",
                      "pacman", "gargle", "ids", "systemfonts", "textshaping",
                      "googledrive", "googlesheets4", "httr", "ragg", "rvest",
                      "xml2", "textshaping", "ragg", "gtools")

new.packages <- list.of.packages[!(list.of.packages %in%
installed.packages()[,"Package"])]

if (length(new.packages) > 0) {
  install.packages(new.packages,
                   dependencies = TRUE,
                   repos = "https://cran.r-project.org")
}

list.of.BiocPackages <- c("DESeq2", "BiocParallel")
new.BiocPackages <- list.of.BiocPackages[!(list.of.BiocPackages %in%
                                             installed.packages()[,"Package"])]

if (length(new.BiocPackages) > 0) {
  BiocManager::install(new.BiocPackages, dependencies = TRUE, ask = FALSE)
}

pacman::p_load(optparse, tidyverse, stringr, DESeq2, tibble, BiocParallel,
               furrr, future, future.batchtools, vsn)