########################
### 2023-07-31 #########
### H. Polveche ########
########################


###################
### trash tests ###
###################


library(tidyverse)

exons <- read.csv2("../exons_genomiques_bis.csv", header = T)
last_exons <- exons %>% 
  group_by(id_gene) %>% 
  dplyr::mutate(
    first = dplyr::first(pos_sur_gene),
    last = dplyr::last(pos_sur_gene)
  ) %>% 
  filter(pos_sur_gene == last) %>% 
  select(id_gene, last, longueur) %>% 
  ungroup() 


inf50 <- last_exons %>% 
  filter(longueur < 50 )

inf100 <- last_exons %>% 
  filter(longueur < 100 )

resu <- last_exons %>% 
  summarise(mean=mean(longueur), sd = sd(longueur), 
            min = min(longueur), max= max(longueur))

# commun_expression <- function(words){
#   ### words : string vector 
#   ### find common expression
#   
#   words.split <- strsplit(words, '')
#   words.split <- lapply(words.split, `length<-`, max(nchar(words)))
#   words.mat <- do.call(rbind, words.split)
#   common.substr.length <- which.max(apply(words.mat, 2, function(col) !length(unique(col)) == 1)) - 1
#   reg.ex <- substr(words[1], 1, common.substr.length)
#   
#   return(reg.ex) 
# }

# c("NCKAP5L_ENSG00000167566", "SH3TC1_ENSG00000125089", "NCS1_ENSG00000107130",
#   "RAB36_ENSG00000100228", "PPARD_ENSG00000112033", "FBLN1_ENSG00000077942")


library(tidyverse)

DEG.ctrl <-  read.csv2("./results/siGL2/normcounts.csv", header = T)
DEG.ctrl <- data.frame(ID=NA, DEG.ctrl)
DEG.ctrl$ID <- rownames(DEG.ctrl)

coldata.ctrl <- read.csv2("./results/siGL2/coldata.csv", header = T)
rownames(coldata.ctrl) <- str_replace(rownames(coldata.ctrl), "-", ".")

coldata.siDDX <- read.csv2("./results/siDDX/coldata.csv", header = T)
coldata.siDDX$sampleName <- str_replace(coldata.siDDX$sampleName, "-", ".")

DEG.siDDX <- read.csv2("./results/siDDX/normcounts.csv", header = T)
DEG.siDDX <- data.frame(ID=NA, DEG.siDDX)
DEG.siDDX$ID <- rownames(DEG.siDDX)

DEG.siDDX.sum <- DEG.siDDX %>% 
  pivot_longer(!ID, names_to = "samples", values_to = "count") %>% 
  mutate(condition = if_else(str_ends(samples, "_RT"), "RT", "canonical")) %>% 
  group_by(condition, ID) %>% 
  summarise(mean = mean(count)) %>% 
  pivot_wider(names_from = c(condition), values_from = mean) 
  
DEG.ctrl.sum <- DEG.ctrl %>% 
  pivot_longer(!ID, names_to = "samples", values_to = "count") %>% 
  mutate(condition = if_else(str_ends(samples, "_RT"), "RT", "canonical")) %>% 
  group_by(condition, ID) %>% 
  summarise(mean = mean(count)) %>% 
  pivot_wider(names_from = c(condition), values_from = mean) 


#DEG.siDDX <- merge(DEG.siDDX, coldata.siDDX, 
#                   by.x = "samples", by.y = "sampleName", all.x = T)



















