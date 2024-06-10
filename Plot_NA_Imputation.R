### Script pour comparer le nombre de Na avant et après l'imputation

setwd("~/Documents/Worms/Plot_GATK/Pop")

library(tidyr)
library(ggplot2)

# Charger les fichiers CSV
files_b <- list.files(pattern = "output_.*_bImputation*")
files_name_b <- sub("output_(.*)_bImputation_PSC.vchk",'\\1',files_b)
dfs_b <- lapply(files_b, read.delim)

files_a <- list.files(pattern = "output_.*_aImputation*")
files_name_a <- sub("output_(.*)_aImputation_PSC.vchk",'\\1',files_a)
dfs_a <- lapply(files_a, read.delim)

Nbr_Na <- data.frame(Population = character(0),
                     NA_b = numeric(0),
                     NA_a = numeric(0))

for(i in 1:length(files)){
  # Set columns name
  colnames(dfs_a[[i]]) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions",	"nIndels","average_depth","nSingleton","nHapRef","nHapAlt","nMissing")
  # Keep only columns we are intersseted in
  dfs_a[[i]] <- dfs_a[[i]][c("sample","nRefHom","nNonRefHom","nHets","nSingleton","nMissing")]
  
  # Set columns name
  colnames(dfs_b[[i]]) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions",	"nIndels","average_depth","nSingleton","nHapRef","nHapAlt","nMissing")
  # Keep only columns we are intersseted in
  dfs_b[[i]] <- dfs_b[[i]][c("sample","nRefHom","nNonRefHom","nHets","nSingleton","nMissing")]
  
  Nbr_Na <- data.frame(rbind(Nbr_Na,data.frame(Population = files_name_a[[i]],
                                         NA_b = mean(dfs_b[[i]]$nMissing),
                                         NA_a = mean(dfs_a[[i]]$nMissing))))
}

library(ggplot2)

# Calculer la différence entre les valeurs "before imputation" et "after imputation"
Nbr_Na$Difference <- Nbr_Na$NA_b - Nbr_Na$NA_a

# Créer le graphique avec les deux distributions et la différence affichée au-dessus des points rouges
ggplot(data = Nbr_Na, aes(x = Population)) +
  geom_point(aes(y = NA_b, color = 'Before Imputation')) +
  geom_point(aes(y = NA_a, color = 'After Imputation')) +
  geom_text(data = subset(Nbr_Na, !is.na(Difference)), aes(y = NA_a, label = paste('+',round(Difference))), vjust = -1.5) +
  labs(y = "log(NA Values)", x = "Population", color = "Imputation") +
  scale_color_manual(values = c('Before Imputation' = 'blue', 'After Imputation' = 'red')) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_y_log10()  # Utiliser une échelle logarithmique 

ggsave('Impact_Imputation.pdf')
