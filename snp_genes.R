library(biomaRt)
library(dplyr)
library(tidyr)
library(gridExtra)

# setwd("/Users/reikotachibana/Documents/Chung Lab/Chimerism") 
# setwd("chimerism_analysis")

# Get genes and their chromosomes/positions
# ensembl <- useMart("ensembl")
ensembl <- useMart("ensembl", host = "https://useast.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# listAttributes(ensembl)
genes <- getBM(attributes = c("ensembl_gene_id", 
                              "external_gene_name",
                              "chromosome_name",
                              "start_position", 
                              "end_position"),
               mart = ensembl)

# Find which gene each SNP is
bed_files <- list.files("DonorHostSNPs", full.names = TRUE)
# ppid_folders <- list.files(pattern = "PPID")
agg_df <- data.frame()
for (bed_file_path in bed_files){
  # bed_file_path <- paste0(ppid, "/", ppid, "_Allvars.bed") # edit
  if (file.exists(bed_file_path)) {
    bed <- read.table(bed_file_path)
    bed$V1 <- sub("^chr", "", bed$V1)
    bed_gr <- GRanges(seqnames = bed$V1, 
                      ranges = IRanges(start = bed$V2, end = bed$V3))
    
    gene_gr <- GRanges(seqnames = genes$chromosome_name,
                       ranges = IRanges(start = genes$start_position,
                                        end = genes$end_position))
    overlaps <- findOverlaps(bed_gr, gene_gr)
    overlapping_genes <- genes[subjectHits(overlaps), ]
    
    gene_names <- overlapping_genes$external_gene_name
    bed$gene <- NA
    bed$gene[queryHits(overlaps)] <- gene_names
    
    file_name <- basename(bed_file_path)
    ppid <- sub("_DonorHostSNPs.bed", "", file_name)
    bed$PPID <- ppid
    
    agg_df <- rbind(agg_df, bed)
    # print(ppid)
    # print(table(bed$gene))
  }
}

plots_list <- list()
snp_types <- c("HET.WT", "HOM.WT", "HET.HOM", "HOM.HET", "WT.HOM")


# HLA
agg_df <- agg_df[grepl("HLA", agg_df$gene), ]
for (snp_type in snp_types){
  filtered_agg_df <- agg_df[agg_df$V5 == snp_type, ]
  agg_summary <- filtered_agg_df %>%
    group_by(PPID, gene) %>%
    summarise(SNP_count = n(), .groups = 'drop')
  agg_wide <- agg_summary %>%
    pivot_wider(names_from = PPID, values_from = SNP_count, values_fill = 0)
  # agg_normalized <- agg_wide %>%
  #   mutate(across(-gene, ~ . / sum(.)))
  agg_long <- agg_wide %>%
    pivot_longer(cols = -gene, names_to = "PPID", values_to = "SNP_count")
  
  plot <- ggplot(agg_long, aes(x = gene, y = PPID, fill = SNP_count)) +  
    geom_tile(color = "white") + 
    scale_fill_gradient(low = "lightblue", high = "darkblue") + 
    labs(title = snp_type, 
         x = "Gene",
         y = "PPID",  
         fill = "SNP Frequency") +  
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6, face = "bold"))  
  plots_list[[snp_type]] <- plot
}
do.call(grid.arrange, c(plots_list, ncol = 2))


for (snp_type in snp_types){
  filtered_agg_df <- agg_df[agg_df$V5 == snp_type, ]
  agg_summary <- filtered_agg_df %>%
    group_by(PPID, gene) %>%
    summarise(SNP_count = n(), .groups = 'drop')
  agg_wide <- agg_summary %>%
    pivot_wider(names_from = PPID, values_from = SNP_count, values_fill = 0)
  # agg_normalized <- agg_wide %>%
  #   mutate(across(-gene, ~ . / sum(.)))
  agg_long <- agg_wide %>%
    pivot_longer(cols = -gene, names_to = "PPID", values_to = "SNP_count")
  
  if (snp_type == "WT.HOM"){
    agg_filtered <- agg_long %>%
      group_by(gene) %>%
      filter(!all(SNP_count <= 2)) %>%  # Keep genes that do not have all counts equal to 1
      ungroup()
  }else if (snp_type == "HOM.WT" | snp_type == "HOM.HET"){
    agg_filtered <- agg_long %>%
      group_by(gene) %>%
      filter(!all(SNP_count <= 5)) %>%  # Keep genes that do not have all counts equal to 1
      ungroup()
  } else if (snp_type == "HET.HOM"){
    agg_filtered <- agg_long %>%
      group_by(gene) %>%
      filter(!all(SNP_count <= 5)) %>%  # Keep genes that do not have all counts equal to 1
      ungroup()
  } else{
    agg_filtered <- agg_long %>%
      group_by(gene) %>%
      filter(!all(SNP_count <= 10)) %>%  # Keep genes that do not have all counts equal to 1
      ungroup()
  }
  
  plot <- ggplot(agg_filtered, aes(x = gene, y = PPID, fill = SNP_count)) +  
            geom_tile(color = "white") + 
            scale_fill_gradient(low = "lightblue", high = "darkblue") + 
            labs(title = snp_type, 
                 x = "Gene",
                 y = "PPID",  
                 fill = "SNP Frequency") +  
            theme_minimal() + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6, face = "bold"))  
  plots_list[[snp_type]] <- plot
}
do.call(grid.arrange, c(plots_list, ncol = 2))










agg_summary <- agg_df %>%
  group_by(PPID, gene) %>%
  summarise(SNP_count = n(), .groups = 'drop')
agg_wide <- agg_summary %>%
  pivot_wider(names_from = PPID, values_from = SNP_count, values_fill = 0)
agg_normalized <- agg_wide %>%
  mutate(across(-gene, ~ . / sum(.)))
agg_long <- agg_normalized %>%
  pivot_longer(cols = -gene, names_to = "PPID", values_to = "SNP_count")
agg_filtered <- agg_long %>%
  group_by(gene) %>%
  filter(!all(SNP_count <= 0.01)) %>%  # Keep genes that do not have all counts equal to 1
  ungroup()

# Create the heatmap with genes on the x-axis and PPIDs on the y-axis
ggplot(agg_filtered, aes(x = gene, y = PPID, fill = SNP_count)) +  # Switch x and y
  geom_tile(color = "white") +  # Add white borders around tiles
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Color gradient
  labs(title = "Heatmap of SNP Counts by Gene and PPID", 
       x = "Gene",  # Change x-axis label to "Gene"
       y = "PPID",  # Change y-axis label to "PPID"
       fill = "SNP Count") +  # Axis and title labels
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))  # Rotate x-axis text for better visibility
