library(MutationalPatterns)
library(BSgenome)
library(vcfR)
library(GenomicRanges)

# setwd("/Users/reikotachibana/Documents/Chung Lab/Chimerism") 

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
# library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)

bed_files <- list.files("DonorHostSNPs", full.names = TRUE)

agg_bed <- data.frame()
for (bed_file_path in bed_files){
  if (file.exists(bed_file_path)) {
    bed <- read.table(bed_file_path)
    bed$REF <- substr(bed$V4, (nchar(bed$V4) - 2), (nchar(bed$V4) - 2))
    bed$ALT <- substr(bed$V4, (nchar(bed$V4)), nchar(bed$V4))
    agg_bed <- rbind(agg_bed, bed)
  }
}


bed_gr <- GRanges(seqnames = agg_bed$V1, 
                  ranges = IRanges(start = agg_bed$V2, width = 1))
bed_gr$REF <- agg_bed$REF
bed_gr$ALT <- agg_bed$ALT

GenomeInfoDb::genome(bed_gr) <- 'hg38'

type_occurrences <- mut_type_occurrences(bed_gr, ref_genome)
# plot_spectrum(type_occurrences, error_bars = 'none')

context <- mut_context(bed_gr, ref_genome)

mut_mat <- mut_matrix(vcf_list = bed_gr, ref_genome = ref_genome)
plot_96 <- plot_96_profile(mut_mat, ymax = 0.055, condensed = TRUE)
ggsave("Images/MutationalPatterns.png", plot = plot_96, width = 16, height = 9)

norm_mut_matrix <- apply(mut_mat, 2, function(x) x / sum(x))
# Get substitution and context from rownames and make long.
tb <- norm_mut_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("full_context") %>%
  dplyr::mutate(
    substitution = stringr::str_replace(full_context, "\\w\\[(.*)\\]\\w", "\\1"),
    context = stringr::str_replace(full_context, "\\[.*\\]", "\\.")
  ) %>%
  dplyr::select(-full_context) %>%
  tidyr::pivot_longer(c(-substitution, -context), names_to = "sample", values_to = "freq") %>% 
  dplyr::mutate(sample = factor(sample, levels = unique(sample)))

colors <- c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE"
)
mutations_plot <- ggplot(data = tb, aes(
                          x = context,
                          y = freq,
                          fill = substitution,
                          width = 0.6
                        )) +
                          geom_bar(stat = "identity", colour = "black", size = .2) +
                          scale_fill_manual(values = colors) +
                          facet_grid(sample ~ substitution) +
                          ylab("Relative contribution") +
                          scale_y_continuous(limits = c(0, 0.055), 
                                             expand = c(0, 0), breaks = NULL) + 
                          guides(fill = "none") +
                          theme_bw() +
                          theme(
                            axis.text.y = element_blank(),  
                            axis.title.x = element_text(size = 12),
                            axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5),
                            panel.grid.major.x = element_blank(),
                            # panel.spacing.x = unit(0.5, "lines"),
                            strip.text.y = element_blank()
                          )
ggsave("Images/MutationalPatterns_2.png", plot = mutations_plot, width = 16, height = 9)
