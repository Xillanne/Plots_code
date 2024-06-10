setwd("~/Documents/Worms/Plot_GATK/Pop")


# Charger les fichiers vcf
files <- list.files(pattern = ".*_final.vcf.gz")

library(vcfR)

dfs <- lapply(files, read.vcfR)
files_name <- sub('(.*)_final.vcf.gz',"\\1",files)



##### First step : kinship matrix ######

for(i in dfs){
  vcf<-dfs[[i]]
  
  # Extract genotype information
  genotype_info <- vcf@gt
  
  # Extract SNP positions
  snp_positions <- vcf@fix[, c("CHROM", "POS")]
  
  write.csv(snp_positions,paste0(files_name,'_snp_positions.csv'),row.names = FALSE,quote = FALSE)
}

