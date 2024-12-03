  ppid <- "PPID1160"

library(data.table)
library(vcfR)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readxl)

#######################
## Donor Chimerism ###
#######################
# setwd("/Users/reikotachibana/Documents/Chung Lab/chimerism_analysis") 
# setwd("C:/Users/Reiko/Documents/chimerism_analysis")
source('chimerism_functions.R')

metadata <- read_excel("ChimerismFileMetaData.xlsx")

pre_fileName <- metadata %>%
  filter(PPID == ppid, (Population == "BULK" | Population == "90- PROG"), `Tx Status` == "Pre") %>%
  pull(`File Name`)
preSNP_file <- paste0(
                      ppid, 
                      "/",
                      pre_fileName,
                      ".snps_filtered.vcf.gz")
preVCF_file <- paste0(ppid, 
                      "/",
                      pre_fileName,
                      ".vcf.gz")

if (grepl("SC", ppid)){
 post_fileName <- metadata %>%
   filter(PPID == ppid, Population == "BULK", `Tx Status` == "Post") %>%
   pull(`File Name`)
} else{
  metadata <- metadata %>%
    mutate(HMTB_Num = as.numeric(gsub("HMTB", "", HMTB)))
  
  post_fileName <- metadata %>%
    filter(PPID == ppid, Population == "BULK", `Tx Status` == "Post") %>%
    slice_min(HMTB_Num) %>%
    pull(`File Name`)
}

postSNP_file <- paste0(ppid, 
                       "/",
                       post_fileName,
                       ".snps_filtered.vcf.gz")
postVCF_file <- paste0(ppid, 
                       "/",
                       post_fileName,
                       ".vcf.gz")

pre_fileName
post_fileName

#################
## WRITE Vars ###
#################
allVafs.cov <- find_allVafs(preSNP_file, preVCF_file, postSNP_file, postVCF_file)
write.table(allVafs.cov, paste0("DonorHostSNPs/", ppid, "_DonorHostSNPs.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
table(allVafs.cov$DonorHost)