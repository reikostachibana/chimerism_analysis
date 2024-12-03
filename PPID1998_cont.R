library(ggplot2)
library(readxl)
library(gridExtra)

setwd("/Users/reikotachibana/Documents/Chung Lab/Chimerism/PPID1998")
source("../chimerism_functions.R")

# Identify HOM shared by both host and donor 

vcf_preBulk <- process_VCF("SNPs/1799FL-23-01-93_S93_L001.snps_filtered.vcf.gz")
HOM_preBulk <- filter_vcf(vcf_preBulk, 30, 1, 1, TRUE)

vcf_postBulk <- process_VCF("SNPs/1799FL-21-01-94_S94_L003.snps_filtered.vcf.gz") # dim==10739
HOM_postBulk <- filter_vcf(vcf_postBulk, 30, 0.9, 1, TRUE)

sharedSNP_postBulk <- semi_join(HOM_postBulk, HOM_preBulk,
                                by = c("CHROM", "POS", "ALT"))
sharedSNP_postBulk <- find_vaf(sharedSNP_postBulk)

# Density plot of VAFs in each populations

# ggplot(sharedSNP_postBulk, aes(x = VAF, y = after_stat(scaled))) +
#   geom_density()

metadata <- read_xlsx("PPID1998.xlsx")
metadata <- as.data.frame(metadata)
unique_hmtb <- unique(metadata$HMTB)
hmtb_list <- as.list(unique_hmtb)

population_colors <- c(
  "ELSE" = "purple", 
  "BULK" = "lightgreen", 
  "HSC" = "darkgreen", 
  "LMPP" = "lightblue", 
  "PROG" = "darkblue"
)

plots_list <- list()
for (hmtb in hmtb_list){
  metadata_filt <- metadata[metadata$HMTB == hmtb, ]
  
  populationDf_list <- list()
  for (i in 1:nrow(metadata_filt)){
    row <- metadata_filt[i, ]
    population <- row$Population
    file_path <- paste0("SNPs/", row$"File Name")
    vcf <- process_VCF(file_path)
    SNPs <- semi_join(vcf, sharedSNP_postBulk,
                      by = c("CHROM", "POS", "ALT"))
    SNPs <- find_vaf(SNPs)
    
    populationDf_list[[population]] <- SNPs
  }
  
  plot <- ggplot(sharedSNP_postBulk, aes(x = VAF, y = after_stat(scaled))) + 
    labs(title = hmtb)
  for (population in names(populationDf_list)) {
    plot <- plot + geom_density(data = populationDf_list[[population]], 
                                aes(x = VAF), fill = population_colors[[population]], 
                                alpha = 0.5)
  }
  plot <- plot +
    theme_minimal()
  plots_list[[hmtb]] <- plot
}
do.call(grid.arrange, c(plots_list, ncol = 2))




vcf_postHSC <- process_VCF("SNPs/1799FL-21-01-89_S89_L003.snps_filtered.vcf.gz")
SNP_postHSC <- semi_join(vcf_postHSC, sharedSNP_postBulk,
                         by = c("CHROM", "POS", "ALT"))
SNP_postHSC <- find_vaf(SNP_postHSC)

vcf_postLMPP <- process_VCF("SNPs/1799FL-21-01-88_S88_L003.snps_filtered.vcf.gz")
SNP_postLMPP <- semi_join(vcf_postLMPP, sharedSNP_postBulk,
                         by = c("CHROM", "POS", "ALT"))
SNP_postLMPP <- find_vaf(SNP_postLMPP)

vcf_postPROG <- process_VCF("SNPs/1799FL-21-01-93_S93_L003.snps_filtered.vcf.gz")
SNP_postPROG <- semi_join(vcf_postPROG, sharedSNP_postBulk,
                          by = c("CHROM", "POS", "ALT"))
SNP_postPROG <- find_vaf(SNP_postPROG)

combined_data <- bind_rows(sharedSNP_postBulk, SNP_postHSC, SNP_postLMPP, SNP_postPROG)

ggplot(sharedSNP_postBulk, aes(x = VAF, y = after_stat(scaled))) +
  geom_density(data = sharedSNP_postBulk, aes(x = VAF), fill = "green", alpha = 0.5) +
  geom_density(data = SNP_postHSC, aes(x = VAF), fill = "blue", alpha = 0.5) +
  geom_density(data = SNP_postLMPP, aes(x = VAF), fill = "lightblue", alpha = 0.5) +
  geom_density(data = SNP_postPROG, aes(x = VAF), fill = "red", alpha = 0.5)

