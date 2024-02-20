if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("Biostrings", "ShortRead", "doParallel",
                       "GenomicAlignments", "Gviz", "GenomicFeatures",
                       "Rsubread"))