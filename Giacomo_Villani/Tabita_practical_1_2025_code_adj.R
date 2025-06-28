#remember to set wd

## reading fst of chr 2
AFR_EUR<- read.table("AFR_EUR.weir.fst", header = T)
AFR_EAS<- read.table("AFR_EAS.weir.fst", header = T) 
EAS_EUR<- read.table("EAS_EUR.weir.fst", header = T)

## no duplicates
AFR_EUR_unique <- AFR_EUR[!duplicated(AFR_EUR$POS), ]
AFR_EAS_unique <- AFR_EAS[!duplicated(AFR_EAS$POS), ]
EAS_EUR_unique <- EAS_EUR[!duplicated(EAS_EUR$POS), ]

## no missing data
AFR_EUR_unique <- AFR_EUR_unique[which(!is.na(AFR_EUR_unique$WEIR_AND_COCKERHAM_FST)), ]
AFR_EAS_unique <- AFR_EAS_unique[which(!is.na(AFR_EAS_unique$WEIR_AND_COCKERHAM_FST)), ]
EAS_EUR_unique <- EAS_EUR_unique[which(!is.na(EAS_EUR_unique$WEIR_AND_COCKERHAM_FST)), ]

## negatives to 0
AFR_EUR_unique$WEIR_AND_COCKERHAM_FST[AFR_EUR_unique$WEIR_AND_COCKERHAM_FST < 0] <- 0
AFR_EAS_unique$WEIR_AND_COCKERHAM_FST[AFR_EAS_unique$WEIR_AND_COCKERHAM_FST < 0] <- 0
EAS_EUR_unique$WEIR_AND_COCKERHAM_FST[EAS_EUR_unique$WEIR_AND_COCKERHAM_FST < 0] <- 0

## The datasets now have different sets of SNPs, filter the files by overlapping the SNPs between them.
shared_pos <- Reduce(intersect, list(
  AFR_EUR_unique$POS,
  AFR_EAS_unique$POS,
  EAS_EUR_unique$POS
))

AFR_EUR_filtered <- AFR_EUR_unique[AFR_EUR_unique$POS %in% shared_pos, ]
AFR_EAS_filtered <- AFR_EAS_unique[AFR_EAS_unique$POS %in% shared_pos, ]
EAS_EUR_filtered <- EAS_EUR_unique[EAS_EUR_unique$POS %in% shared_pos, ]


POS_target <- 109513601

AFR_EUR_filtered[AFR_EUR_filtered$POS==POS_target,]
AFR_EAS_filtered[AFR_EAS_filtered$POS==POS_target,]
EAS_EUR_filtered[EAS_EUR_filtered$POS==POS_target,]

AFR_EUR_filtered <- AFR_EUR_filtered[order(AFR_EUR_filtered$WEIR_AND_COCKERHAM_FST), ]
AFR_EAS_filtered <- AFR_EAS_filtered[order(AFR_EAS_filtered$WEIR_AND_COCKERHAM_FST), ]
EAS_EUR_filtered <- EAS_EUR_filtered[order(EAS_EUR_filtered$WEIR_AND_COCKERHAM_FST), ]


## computing PBS based on Fst
Fst_BC <- AFR_EUR_filtered[AFR_EUR_filtered$POS==POS_target,]$WEIR_AND_COCKERHAM_FST
Fst_AC <- AFR_EAS_filtered[AFR_EAS_filtered$POS==POS_target,]$WEIR_AND_COCKERHAM_FST
Fst_AB <- EAS_EUR_filtered[EAS_EUR_filtered$POS==POS_target,]$WEIR_AND_COCKERHAM_FST

PBS_col <- ( (-log(1-Fst_AB)) + (-log(1-Fst_AC)) - (-log(1-Fst_BC)) ) / 2
POS_AE <- AFR_EUR_filtered$POS
PBS <- cbind(POS_AE, PBS_col)
colnames(PBS) <- c("POS", "PBS")


######################## second practical ########################
install.packages("rehh")
library(rehh)
library(tidyverse)

AFR_ehh <- read.table("Chr2_EDAR_LWK_500K.recode.vcf") # (African population)
EAS_ehh <- read.table("Chr2_EDAR_CHS_500K.recode.vcf") # (East Asian population)

snp_of_interest <- "rs3827760"
AFR_ehh <- data2haplohh(hap_file="Chr2_EDAR_LWK_500K.recode.vcf", polarize_vcf = F)
EAS_ehh <- data2haplohh(hap_file="Chr2_EDAR_CHS_500K.recode.vcf", polarize_vcf = F)

AFR_ehh_calc <- calc_ehh(haplohh=AFR_ehh, mrk=snp_of_interest)
EAS_ehh_calc <- calc_ehh(haplohh=EAS_ehh, mrk=snp_of_interest)

ggplot() +
  geom_line(data = AFR_ehh_calc[[3]], aes(x = POSITION, y = EHH_A), colour = "red") +
  geom_line(data = EAS_ehh_calc[[3]], aes(x = POSITION, y = EHH_A), colour = "blue") +
  theme_minimal()

furcation_afr_snp1 <- calc_furcation(haplohh = AFR_ehh, mrk = snp_of_interest)
plot(furcation_afr_snp1)
furcation_eas_snp1 <- calc_furcation(haplohh = EAS_ehh, mrk = snp_of_interest)
plot(furcation_eas_snp1)

snp_of_interest <- "rs190120314"
AFR_ehh_calc <- calc_ehh(haplohh=AFR_ehh, mrk=snp_of_interest)
EAS_ehh_calc <- calc_ehh(haplohh=EAS_ehh, mrk=snp_of_interest)

ggplot() +
  geom_line(data = AFR_ehh_calc[[3]], aes(x = POSITION, y = EHH_A), colour = "red") +
  geom_line(data = EAS_ehh_calc[[3]], aes(x = POSITION, y = EHH_A), colour = "blue") +
  theme_minimal()

furcation_afr_snp1 <- calc_furcation(haplohh = AFR_ehh, mrk = snp_of_interest)
plot(furcation_afr_snp1)
furcation_eas_snp1 <- calc_furcation(haplohh = EAS_ehh, mrk = snp_of_interest)
plot(furcation_eas_snp1)

#iHS is a measure of the amount of extended haplotype homozygosity at a given SNP along the ancestral allele relative to the derived allele. This measure is typically standardized empirically to the distribution of observed iHS scores over a range of SNPs with similar derived allele frequencies.
#calculating for all SNPs
AFR_ehh_scan <- scan_hh(haplohh=AFR_ehh)
EAS_ehh_scan <- scan_hh(haplohh=EAS_ehh)

snp_of_interest <- "rs3827760"

AFR_iHS <- ihh2ihs(scan= AFR_ehh_scan, min_maf = 0.02, freqbin = 0.01)
EAS_iHS <- ihh2ihs(scan= EAS_ehh_scan, min_maf = 0.02, freqbin = 0.01)

ggplot() +
  geom_line(data = AFR_iHS$ihs, aes(x = POSITION, y = IHS), colour = "red") +
  geom_line(data = EAS_iHS$ihs, aes(x = POSITION, y = IHS), colour = "gold")+
  theme_minimal()


## mean over sliding windows
# Create a function to estimate the mean in sliding windows.
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(abs(data[spots[i]:(spots[i] + window - 1)]),na.rm=TRUE)
  }
  return(result)
}

EAS_slide_IHS <- slideFunct(EAS_iHS$ihs$IHS, window = 50, step = 40) 

slidePos <- function(data, window, step){
 total <- length(data)
 spots <- seq(from = 1, to = (total - window + 1), by = step)
 result <- vector(length = length(spots))
 for(i in 1:length(spots)){
 result[i] <- data[spots[i]]
 }
 return(result)
}

EAS_slide_POS <- slidePos(EAS_iHS$ihs$IHS, window = 50, step = 40) 