#!/usr/bin/env Rscript
# list.of.packages <- c("optparse", "dplyr", "tidyr", "stringr", "BiocManager", 
#                       "tidyverse", "tibble", "purrr", "furrr", "future", "vsn",
#                       "pacman")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# 
# if (length(new.packages) > 0) {
#   install.packages(new.packages, dependencies = TRUE, repos = "https://cran.r-project.org")
# }
# 
# list.of.BiocPackages <- c("DESeq2", "BiocParallel")
# new.BiocPackages <- list.of.BiocPackages[!(list.of.BiocPackages %in% installed.packages()[,"Package"])]
# 
# if (length(new.BiocPackages) > 0) {
#   BiocManager::install(new.BiocPackages, dependencies = TRUE, ask = FALSE)
# }

pacman::p_load(optparse, tidyverse, stringr, DESeq2, tibble, BiocParallel,
               furrr, future, future.batchtools, vsn)

# Option parser and loading data:
option_list = list(
  make_option(c("-f", "--fusion"), type="character", default=NULL, 
              help="fusions folder path", metavar="character"),
  make_option(c("-c", "--count"), type="character", default=NULL,
              help="htseq counts folder path.", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL,
              help="path to design table", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=4,
              help="number of threads", metavar="integer"),
  make_option(c("-m", "--memory"), type="integer", default=8,
              help="memory in Gbytes", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Number of threads for DESeq2:
register(MulticoreParam(opt$threads))

# Set the number of threads/workers to use
plan(multisession, workers = opt$threads)

# Set memory limits for parallel workers
options(future.globals.maxSize = ((opt$memory*1000)*1024^2))

# read design input file:
design <- read.table(file = opt$design, header = TRUE, row.names = NULL)
list_sample <- as.list(design$sample)

# Function to parse all file from a list:
parse_fusion <- function(ech) {
  # read fusion file:
  fusion_file <- paste0(ech, "_fusions.tsv")
  df1 <- read.table(file = paste0(opt$fusion, fusion_file))
  colnames(df1) <- c("gene1","gene2","strand1",
                     "strand2","breakpoint1","breakpoint2","site1",
                     "site2","type","split_reads1","split_reads2",
                     "discordant_mates","coverage1","coverage2","confidence",
                     "reading_frame","tags","retained_protein_domains",
                     "closest_genomic_breakpoint1",
                     "closest_genomic_breakpoint2","gene_id1","gene_id2",
                     "transcript_id1","transcript_id2","direction1",
                     "direction2","filters","fusion_transcript",
                     "peptide_sequence","read_identifiers")
  # filter fusion file:
  df1 <- filter(df1, type  == 'deletion/read-through' |
                  type == "deletion/read-through/5'-5'" |
                  type == "deletion/read-through/3'-3'")
  df1 <- df1[df1$strand1 == df1$strand2,]
  df1 <- df1 %>% filter(strand1 %in% c("+/+", "-/-"))
  
  df1 <- df1 %>% select(gene1, gene2, strand1, split_reads1, split_reads2, 
                        discordant_mates, type, gene_id1, gene_id2)
  
  htseq_file <- paste0(ech, ".tsv")
  htseq_df <- read.table(file = paste0(opt$count, htseq_file), header = FALSE)
  colnames(htseq_df) <- c("gene_id", "count")
  
  df1 <- left_join(df1, htseq_df, by = join_by(gene_id1 == gene_id),
                   suffix = c(".x", ".y"))
  colnames(df1)[10] <- "count1"
  df1 <- left_join(df1, htseq_df, by = join_by(gene_id2 == gene_id),
                   suffix = c(".x", ".y"))
  colnames(df1)[11] <- "count2"
  colnames(df1)[3] <- "strand"

  df1 <- df1 %>% dplyr::mutate(ID = paste0(gene1, "_", gene2))

  df2 <- df1 %>% dplyr::group_by(ID) %>%
    dplyr::summarise(reads1=sum(split_reads1),
                     reads2=sum(split_reads2),
                     discordant=sum(discordant_mates)) %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(reads_total = sum(reads1, 
                                  reads2, 
                                  discordant)
           ) %>%
    dplyr::select(ID, reads_total)
  
  colnames(df2)[2] <- ech
  
  df1 <- df1 %>% select(all_of(c("ID", "gene_id1", "gene_id2", 
                                 "split_reads1", "split_reads2", 
                                 "discordant_mates", "type", "count1", "count2",
                                 "strand")))
  
  write.table(df1, file = paste0(opt$output, ech, "_parsed_fusion.csv"), 
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  df2$ID <- as.factor(df2$ID)
  return(df2)
}

# all_results <- map(list_sample, parse_fusion)
all_results <- furrr::future_map(list_sample, parse_fusion)
gc()
concat_res <- purrr::reduce(all_results, full_join, by = "ID")
concat_res[is.na(concat_res)] <- 0

# DESeq2:
## Prepare the dataset:
deseq_df <- as.matrix(concat_res %>% column_to_rownames(var = "ID"))
exp_design <- design %>% column_to_rownames(var = "sample")
exp_design$condition <- as.factor(exp_design$condition)
exp_design$rep <- as.factor((exp_design$rep))

dds <- DESeq2::DESeqDataSetFromMatrix(countData=deseq_df, 
                                      colData=exp_design, 
                                      design=~condition + rep, 
                                      tidy = FALSE)

## Filtering:
keep <- rowSums(counts(dds) > 1) >= 3
dds <- dds[keep,]

dds <- DESeq2::DESeq(dds)

# Plots:
pdf(paste0(opt$output, "rplot.pdf"))
res <- results(dds)
plotMA(res) #, ylim=c(-2,2))
plotDispEsts(dds)

# this gives log2(n + 1)
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))
dev.off()

# save results:
save_deseq <- function(comparaison, dds) {
  tmp <- as.data.frame(DESeq2::results(dds, contrast=c("condition",comparaison)))
  tmp <- tibble::rownames_to_column(tmp, "ID")
  write.table(tmp, 
              file = paste0(opt$output, comparaison[1], "_", comparaison[2], ".tsv"), 
              row.names = FALSE,
              dec = ",",
              sep = "\t",
              qmethod = "double"
  )
}

sapply(list_comparaison, save_deseq, dds)
