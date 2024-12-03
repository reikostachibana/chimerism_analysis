library(GenomicRanges)
library(ggplot2)
library(MutationalPatterns)

# setwd("/Users/reikotachibana/Documents/Chung Lab/Chimerism") 

# ppid_folders <- list.files(pattern = "PPID")
hyb_probes <- read.table("xgen-aml-cancer-hyb-panel-probes_hg38.bed")
hla_probes <- read.table("HLAsiteProbes.bed")
hyb_probes_gr <- GRanges(seqnames = hyb_probes$V1, 
                     ranges = IRanges(start = hyb_probes$V2, end = hyb_probes$V3))
hla_probes_gr <- GRanges(seqnames = hla_probes$V1,
                         ranges = IRanges(start = hla_probes$V2, end = hla_probes$V3))

bed_files <- list.files("DonorHostSNPs", full.names = TRUE)

agg_df <- data.frame()
for (bed_file_path in bed_files){
  if (file.exists(bed_file_path)) {
    bed <- read.table(bed_file_path)
    bed$REF <- substr(bed$V4, (nchar(bed$V4) - 2), (nchar(bed$V4) - 2))
    bed$ALT <- substr(bed$V4, (nchar(bed$V4)), nchar(bed$V4))
    bed_gr <- GRanges(seqnames = bed$V1, 
                      ranges = IRanges(start = bed$V2, end = bed$V3))
    bed_gr$REF <- bed$REF
    bed_gr$ALT <- bed$ALT
    
    hyb_overlaps <- data.frame(findOverlaps(bed_gr, hyb_probes_gr))
    hyb_overlaps$pos <- bed$V2[hyb_overlaps$queryHits]
    hyb_overlaps$start <- hyb_probes$V2[hyb_overlaps$subjectHits]
    hyb_overlaps$relPos <- hyb_overlaps$pos - hyb_overlaps$start + 1
    hyb_overlaps$DonorHost <- bed$V5[hyb_overlaps$queryHits]
    
    hla_overlaps <- data.frame(findOverlaps(bed_gr, hla_probes_gr))
    hla_overlaps$pos <- bed$V2[hla_overlaps$queryHits]
    hla_overlaps$start <- hla_probes$V2[hla_overlaps$subjectHits]
    hla_overlaps$relPos <- hla_overlaps$pos - hla_overlaps$start + 1
    hla_overlaps$DonorHost <- bed$V5[hla_overlaps$queryHits]
    
    combined_overlaps <- rbind(hyb_overlaps, hla_overlaps)
    
    file_name <- basename(bed_file_path)
    ppid <- sub("_DonorHostSNPs.bed", "", file_name)
    combined_overlaps$PPID <- ppid
    
    agg_df <- rbind(agg_df, combined_overlaps)
  } 
}

# type_occurrences <- mut_type_occurrences(grl, ref_genome)

ggplot(agg_df, aes(x = relPos, fill = DonorHost)) +
  ggsci::scale_fill_jama(alpha = .5)+
  geom_density(aes( y = after_stat(density)))+
  ggpubr::theme_pubr()+
  facet_wrap(~PPID, scales = "fixed")
ggsave("Images/probes.png", width = 16, height = 9)


# ggplot(agg_df, aes(x = relPos, fill = PPID)) +
#   ggsci::scale_fill_jama(alpha = .5)+
#   geom_density(aes( y = after_stat(density)))+
#   ggpubr::theme_pubr()+
#   facet_wrap(~DonorHost)
# 
# 
# ggplot(agg_df, aes(x = relPos, fill = DonorHost)) +
#   ggsci::scale_fill_jama(alpha = .5)+
#   geom_density(aes( y = after_stat(density)))+
#   ggpubr::theme_pubr()+
#   facet_wrap(~DonorHost)
