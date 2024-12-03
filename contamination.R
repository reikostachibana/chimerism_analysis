ppid <- "PPID1998"

library(ggplot2)
library(readxl)
library(gridExtra)
library(gtools)

# setwd("/Users/reikotachibana/Documents/Chung Lab/Chimerism/PPID1998")
setwd("chimerism_analysis")
source("chimerism_functions.R")

# Identify HOM shared by both host and donor 

metadata <- read_excel("ChimerismFileMetaData.xlsx")

if (grepl("Farrel", ppid)){
  pre_fileName <- "45_S22"
  preSNP_file <- paste0(ppid, "/", pre_fileName, ".snps_filtered.vcf.gz")
  preVCF_file <- paste0(ppid, "/", pre_fileName, ".vcf.gz")
  
  post_fileName <- "Farrell_30_Bulk_S2"
  postSNP_file <- paste0(ppid, "/", post_fileName, ".snps_filtered.vcf.gz")
  postVCF_file <- paste0(ppid, "/", post_fileName, ".vcf.gz")
} else{
  metadata <- metadata %>%
    mutate(HMTB_Num = as.numeric(gsub("HMTB", "", HMTB)))
  
  pre_fileName <- metadata %>%
    filter(PPID == ppid, (Population == "BULK" | Population == "90- PROG"), `Tx Status` == "Pre") %>%
    slice_min(HMTB_Num) %>%
    pull(`File Name`)
  preSNP_file <- paste0(ppid, "/", pre_fileName, ".snps_filtered.vcf.gz")
  preVCF_file <- paste0(ppid, "/", pre_fileName, ".vcf.gz")
  
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
  
  postSNP_file <- paste0(ppid, "/", post_fileName, ".snps_filtered.vcf.gz")
  postVCF_file <- paste0(ppid, "/", post_fileName, ".vcf.gz")
}

vcf_preBulk <- process_VCF(preSNP_file)
HOM_preBulk <- filter_vcf(vcf_preBulk, 30, 1, 1, TRUE)

vcf_postBulk <- process_VCF(postSNP_file) # dim==10739
HOM_postBulk <- filter_vcf(vcf_postBulk, 30, 0.9, 1, TRUE)

sharedSNP_postBulk <- semi_join(HOM_postBulk, HOM_preBulk,
                                by = c("CHROM", "POS", "ALT"))
sharedSNP_postBulk <- find_vaf(sharedSNP_postBulk)

# Density plot of VAFs in each populations

# ggplot(sharedSNP_postBulk, aes(x = VAF, y = after_stat(scaled))) +
#   geom_density()

metadata <- metadata[metadata$PPID == ppid, ]
metadata <- as.data.frame(metadata)
unique_hmtb <- unique(metadata$HMTB)
hmtb_list <- as.list(unique_hmtb)
hmtb_list <- mixedsort(unlist(hmtb_list))

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
    file_path <- paste0(ppid, "/", row$"File Name", ".snps_filtered.vcf.gz")
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
                                alpha = 0.5, adjust = 0.5)
  }
  plot <- plot +
    theme_minimal()
  plots_list[[hmtb]] <- plot
}

ggsave(paste0("contamination_", ppid, ".png"), plot = do.call(grid.arrange, c(plots_list, ncol = 2)),
       width = 16, height = 9)