### Script to summarize SNPs effects
args <- commandArgs(trailingOnly = TRUE)
file <-  args[1] 
output <- args[2]

library(HDInterval)
library(data.table)
library(ggplot2)
library(matrixStats)


Df <- as.data.frame(fread(file))
print('Data frame loaded')

row.names(Df) <- Df$V1
Df$V1 <- NULL
chrom <- Df$CHROM
Df$CHROM <- NULL
Df$POS <- NULL

HDI<-hdi(t(Df),0.98)
print('Credibility interval computed')

is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}

sig <- colnames(HDI)[!apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"]))]
print('Credible SNPs found')

# Ajouter des colonnes pour la médiane et les intervalles de crédibilité
Df$median <- rowMedians(as.matrix(Df))
Df$lower <- HDI['lower',]
Df$upper <- HDI['upper',]

# Ajouter une colonne pour indiquer les SNPs significatifs
Df$credible <- ifelse(rownames(Df) %in% sig, "Credible", "Not Credible")

Df$CHROM <- as.factor(chrom)

write.csv(Df[,c('CHROM','median','lower','upper','credible')],gsub('XXX',output,'Summary_XXX.csv'))
