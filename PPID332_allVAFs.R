library(data.table)
library(vcfR)
library(ggplot2)
library(dplyr)
library(ggpubr)

#######################
## Donor Chimerism ###
#######################
setwd("/Users/reikotachibana/Documents/Chung Lab") 
source('./Chimerism/chimerism_functions.R')
source('./Chimerism/chimerism_functions_v2.R')


ppid <- "PPID332"
preSNP_file <- paste0(ppid, 
                      "/SNPs/HMTB_680_Bulk_S7.snps_filtered.vcf.gz")
preVCF_file <- paste0(ppid, 
                      "/VCFs/HMTB_680_Bulk_S7.vcf.gz")
postSNP_file <- paste0(ppid,
                       "/SNPs/HMTB_1127_Bulk_S11.snps_filtered.vcf.gz")
postVCF_file <- paste0(ppid, 
                       "/VCFs/HMTB_1127_Bulk_S11.vcf.gz")


#################
## WRITE Vars ###
#################
allVafs.cov <- find_allVafs(preSNP_file, preVCF_file, postSNP_file, postVCF_file)
write.table(allVafs.cov, paste0(ppid, "/", ppid, "_Allvars_v2.bed"), 
            quote = F, sep = "\t", row.names = F, col.names = F)