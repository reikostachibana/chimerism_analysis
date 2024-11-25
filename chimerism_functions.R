library(vcfR)
library(data.table)
library(dplyr)

process_gtdf <- function(gtdf_raw){
  split_values <- strsplit(as.character(gtdf_raw$'3'), ":")
  gtdf <- do.call(rbind, lapply(split_values, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  colnames(gtdf) <- c('GT', 'AD', 'DP', 'GQ', 'PL')
  return(gtdf)
}

process_VCF <- function(file_path){ 
  vcf <- read.vcfR( file_path, verbose = FALSE )
  fxdf <- as.data.frame(getFIX(vcf))
  gtdf_raw <- as.data.frame(vcf@gt)
  gtdf <- process_gtdf(gtdf_raw)
  vcfdf <- cbind(fxdf, gtdf)
  setDT(vcfdf)
  vcfdf <- vcfdf[!ALT %like% ","] # Removes multi-allelic sites
  vcfdf$DP <- as.numeric(vcfdf$DP)
  vcfdf$CHROM_POS <- paste(vcfdf$CHROM,
                           vcfdf$POS,
                           sep = "_")
  vcfdf$BARCODE <- paste(vcfdf$CHROM,
                         vcfdf$POS,
                         vcfdf$REF,
                         vcfdf$ALT,
                         sep = "_")
  return(vcfdf)
}

find_vaf <- function(vcfdf){
  split_AD <- strsplit(as.character(vcfdf$'AD'), ",")
  alt <- as.numeric(sapply(split_AD, function(x) x[2]))
  vcfdf$DP <- as.numeric(vcfdf$DP)
  vcfdf$VAF <- alt / vcfdf$DP
  return(vcfdf)
}

filter_vcf <- function(vcfdf, DP_threshold, VAF_min, VAF_max, apply_filter) {
  if (apply_filter) {
    vcfdf <- vcfdf[vcfdf$FILTER == "PASS", ]
  }
  vcfdf <- vcfdf[vcfdf$DP >= DP_threshold, ]
  vcfdf <- find_vaf(vcfdf)
  vcfdf <- vcfdf[vcfdf$VAF >= VAF_min & vcfdf$VAF <= VAF_max, ]
  return(vcfdf)
}

write_bed <- function(df, filename){
  df$POS <- as.numeric(df$POS)
  df$END <- df$POS + 1
  bed_data <- df[, c("CHROM", "POS", "END")]
  write.table(bed_data, file=filename, sep="\t", quote=FALSE, 
              col.names=FALSE, row.names=FALSE, append=TRUE)
}

read_coverage <- function(filename){
  df <- read.table(filename)
  df <- df %>%
    rename("CHROM" = "V1",
           "POS+1" = "V2",
           "DP" = "V3")
  df$CHROM_POS <- paste(df$CHROM,
                        df$"POS+1" - 1,
                        sep = "_")
  df$POS <- df$"POS+1" - 1
  return (df)
}

find_chimerism <- function(host_SNP, vcf, coverage){ # VCF of post-Tx
  host_SNP_CHROMPOS <- host_SNP$CHROM_POS [(host_SNP$CHROM_POS 
                                            %in% coverage$CHROM_POS)]
  host_SNP_BARCODE <- host_SNP$BARCODE [(host_SNP$CHROM_POS 
                                         %in% coverage$CHROM_POS)]
  vcf_hostSNP <- vcf[vcf$BARCODE %in% host_SNP_BARCODE]
  vcf_hostSNP <- find_vaf(vcf_hostSNP)
  vcf_hostSNP <- vcf_hostSNP[, c("CHROM", "POS", "VAF")]
  vcf_hostSNP$POS <- as.numeric(vcf_hostSNP$POS)
  coverage$POS <- as.numeric(coverage$POS)
  coverage <- coverage[coverage$CHROM_POS %in% host_SNP$CHROM_POS,]
  vcf <- merge(vcf_hostSNP, coverage, by= c("CHROM", "POS"), all.y = T)
  vcf$VAF[is.na(vcf$VAF)] <- 0 
  sum_alt_depth <- sum(vcf$VAF * vcf$DP )
  sum_DP <- sum(vcf$DP)
  return (c((sum_alt_depth / sum_DP), sum_DP))
}

find_all_chimerism <- function(samples, coverages){
  final.data <- data.frame()
  for (i in seq_along(samples)) {
    sample <- samples[i]
    coverage_file <- coverages[i]
    vcf <- process_VCF(sample)
    coverage <- read_coverage(coverage_file)
    c_chimerism <- find_chimerism(host_SNP, vcf, coverage)
    chimerism <- c_chimerism[1]
    dp <- c_chimerism[2]
    df.temp <- data.frame(sample, chimerism, dp)
    final.data <- rbind(final.data, df.temp)
  }
  return()
}


find_allVafs <- function(preSNP_file, preVCF_file, postSNP_file, postVCF_file){
  pre_SNP <- process_VCF(preSNP_file) # dim=
  post_SNP <- process_VCF(postSNP_file) # dim=
  pre_SNP.hom <-filter_vcf(pre_SNP, 50, 1, 1, TRUE) # dim=
  pre_SNP.het <- filter_vcf(pre_SNP, 50, 0.4, 0.6, TRUE)
  post_SNP <-filter_vcf(post_SNP, 50, 0.2, 1, TRUE) # dim=
  
  pre_VCF <- process_VCF(preVCF_file) # dim=
  post_VCF <- process_VCF(postVCF_file) # dim=
  post_VCF.vaf <- find_vaf(post_VCF)
  
  pre_VCF <- filter_vcf(pre_VCF, 0, 0, 1, FALSE)
  pre_VCFhom <- pre_VCF[pre_VCF$VAF > .9,]
  pre_VCFhet <- pre_VCF[(pre_VCF$VAF > .35 & pre_VCF$VAF < .65),]
  post_VCFhom <- filter_vcf(post_VCF, 0, 0.25, 1, FALSE)
  post_VCFhet <- filter_vcf(post_VCF, 0, 0.10, 1, FALSE)
  
  host.hom <- pre_SNP.hom[!(pre_SNP.hom$BARCODE %in% post_VCFhom$BARCODE),]
  host.het <- pre_SNP.het[!(pre_SNP.het$BARCODE %in% post_VCFhet$BARCODE),]
  
  post_SNPsdonor <- post_SNP[!(post_SNP$BARCODE %in% pre_VCF$BARCODE),]
  post_SNPs.HostHom <- post_SNP[(post_SNP$BARCODE %in% pre_VCFhom$BARCODE),]
  post_SNPs.HostHet <- post_SNP[(post_SNP$BARCODE %in% pre_VCFhet$BARCODE),]
  
  post_VCF.host.vafs.hom <- post_VCF.vaf[post_VCF.vaf$BARCODE %in% host.hom$BARCODE,]
  post_VCF.host.vafs.het <- post_VCF.vaf[post_VCF.vaf$BARCODE %in% host.het$BARCODE,]
  
  p1 <- ggplot(post_SNPsdonor, aes(x=VAF, y = after_stat(scaled))) +
    geom_density(data = post_SNPs.HostHom, aes(x=VAF),  fill = "darkgreen", size = 1)+
    geom_density(data = post_SNPs.HostHet, aes(x=VAF),  fill = "lightgreen", size = 1)+
    geom_density(fill = "darkred", size = 1)+
    geom_density(data = post_VCF.host.vafs.hom, aes(x=VAF),  fill = "darkblue", size = 1)+
    geom_density(data = post_VCF.host.vafs.het, aes(x=VAF),  fill = "lightblue", size = 1)+
    # geom_vline(xintercept = .5, linetype = "dashed", size = 2, )+
    # geom_vline(xintercept = 1, linetype = "dashed", size =2 )+
    xlim(c(0,1.1))+
    theme_pubr()
  
  modevafs <- data.frame()
  for (i in 1:3){
    dens <- layer_data(p1, i)
    dens <- dens[order(dens$x),]
    # Run length encode the sign of difference
    rle <- rle(diff(as.vector(dens$y)) > 0)
    # Calculate startpoints of runs
    starts <- cumsum(rle$lengths) - rle$lengths + 1
    # Take the points where the rle is FALSE (so difference goes from positive to negative) 
    maxima_id <- starts[!rle$values]
    maxima <- dens[maxima_id,]
    maxima <- maxima[order(-maxima$y),]
    maxima <- maxima[1:2,]
    maxima <- maxima[order(-maxima$x),]
    maxima$DonorSNP <- "HOM"
    maxima$DonorSNP[2] <- "HET"
    modevafs <- rbind(modevafs, maxima)
    
  }
  modevafs$HostSNP <- "WT"
  modevafs$HostSNP[modevafs$fill == "lightgreen"] <- "HET"
  modevafs$HostSNP[modevafs$fill == "darkgreen"] <- "HOM"
  modevafs$DonorHost <- paste0(modevafs$DonorSNP, ".", modevafs$HostSNP)
  
  modevafs$Vafmin <- modevafs$x-.05
  modevafs$Vafmax <- modevafs$x+.05
  
  post_SNPsdonor$TYPE <- "DonorSpecific"
  post_SNPs.HostHet$TYPE <- "DonorHostHet"
  post_SNPs.HostHom$TYPE <- "DonorHostHom"
  host.hom$TYPE <- "HostHom"
  host.het$TYPE <- "HostHet"
  
  post_SNPsdonor$HostSNP <- "WT" 
  post_SNPsdonor$DonorSNP <- "None"
  post_SNPsdonor$DonorSNP[(post_SNPsdonor$VAF > modevafs$Vafmin[modevafs$DonorHost == "HET.WT"] &
                             post_SNPsdonor$VAF < modevafs$Vafmax[modevafs$DonorHost == "HET.WT"])] <- "HET"
  post_SNPsdonor$DonorSNP[(post_SNPsdonor$VAF > modevafs$Vafmin[modevafs$DonorHost == "HOM.WT"] &
                             post_SNPsdonor$VAF < modevafs$Vafmax[modevafs$DonorHost == "HOM.WT"])] <- "HOM"
  post_SNPs.HostHet$HostSNP <- "HET"
  post_SNPs.HostHet$DonorSNP <- "None"
  post_SNPs.HostHet$DonorSNP[(post_SNPs.HostHet$VAF > modevafs$Vafmin[modevafs$DonorHost == "HET.HET"] &
                                post_SNPs.HostHet$VAF < modevafs$Vafmax[modevafs$DonorHost == "HET.HET"])] <- "HET"
  post_SNPs.HostHet$DonorSNP[(post_SNPs.HostHet$VAF > modevafs$Vafmin[modevafs$DonorHost == "HOM.HET"] &
                                post_SNPs.HostHet$VAF < modevafs$Vafmax[modevafs$DonorHost == "HOM.HET"])] <- "HOM"
  post_SNPs.HostHom$HostSNP <- "HOM"
  post_SNPs.HostHom$DonorSNP <- "None"
  post_SNPs.HostHom$DonorSNP[(post_SNPs.HostHom$VAF > modevafs$Vafmin[modevafs$DonorHost == "HET.HOM"] &
                                post_SNPs.HostHom$VAF < modevafs$Vafmax[modevafs$DonorHost == "HET.HOM"])] <- "HET"
  post_SNPs.HostHom$DonorSNP[(post_SNPs.HostHom$VAF > modevafs$Vafmin[modevafs$DonorHost == "HOM.HOM"] &
                                post_SNPs.HostHom$VAF < modevafs$Vafmax[modevafs$DonorHost == "HOM.HOM"])] <- "HOM"
  
  
  host.hom$HostSNP <- "HOM"
  host.hom$DonorSNP <- "WT"
  host.het$HostSNP <- "HET"
  host.het$DonorSNP <- "WT"
  
  allVafs <- rbind(post_SNPsdonor, post_SNPs.HostHom, post_SNPs.HostHet, host.hom, host.het)
  head(allVafs)
  allVafs$DonorHost <- paste0(allVafs$DonorSNP, ".", allVafs$HostSNP)
  allVafs <- allVafs[!(allVafs$DonorHost %in% c("HET.HET", "HOM.HOM")),]
  allVafs <- allVafs[!(allVafs$DonorHost %like% c("None")),]
  
  
  allVafs.cov <- allVafs[,c("CHROM", "POS", "BARCODE", "DonorHost")]
  allVafs.cov$POS <- as.numeric(allVafs.cov$POS)
  allVafs.cov$POS2 <- allVafs.cov$POS+1
  allVafs.cov <- allVafs.cov[,c("CHROM", "POS", "POS2", "BARCODE", "DonorHost")]
  
  return(allVafs.cov)
}