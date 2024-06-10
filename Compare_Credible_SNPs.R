### Script which compare the credible SNP for each trait
output <- 'VanRaden_A6_NGM_0.99'
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned")
#setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")

library(tidyr)
library(ggplot2)

# Charger les fichiers CSV
files <- list.files(pattern = "Summary_.*_VanRaden_A6_NGM_0.99.csv")
file_names <- sub("Summary_(.*)_VanRaden_A6_NGM_0.99.csv", "\\1", files)
dfs <- lapply(files, read.csv)
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")

# Filtrer les lignes "Credible"
credibles <- lapply(dfs, function(df) df[df$credible == "Credible",]$X)

# Fonction pour comparer les noms de lignes crédibles entre des groupes de dataframes
compare_lists <- function(lists) {
  if (length(lists) == 1) {
    return(lists[[1]])
  } else {
    result <- lists[[1]]
    for (i in 2:length(lists)) {
      result <- intersect(result, lists[[i]])
    }
    return(result)
  }
}



# Générer toutes les combinaisons possibles
library(combinat)

results_solo <- list()
results_solo[[paste0("combinations_", 1)]] <- lapply(seq_along(dfs), function(idx) {
  list(
    combination = file_names[idx],
    intersection = length(credibles[[idx]])
  )
})

results <- list()
for (i in 2:length(credibles)) {
  comb <- combn(seq_along(credibles), i, simplify = FALSE)
  comb_names <- combn(file_names , i, simplify = FALSE)
  results[[paste0("combinations_", i)]] <- lapply(seq_along(comb), function(idx) {
    list(
      combination = comb_names[[idx]],
      intersection = compare_lists(credibles[comb[[idx]]])
    )
  })
}


# Affichage des résultats
# Ouvrir un fichier en écriture
output_file <- file(gsub('XXX',output,'Shared_credible_SNPs_XXX.tsv'), open = "wt")

# Écrire l'entête dans le fichier
cat("Comparison\tNbr_SNPs\n", file = output_file)

# Parcourir les résultats et écrire dans le fichier
for (i in seq_along(results_solo)) {
  for (k in seq_along(results_solo[[i]])) {
    cat(paste(results_solo[[i]][[k]]$combination, collapse = ","), "\t", file = output_file)
    cat(results_solo[[i]][[k]]$intersection, "\n", file = output_file)
  }
  cat("\n", file = output_file)
}

# Parcourir les résultats et écrire dans le fichier
# Parcourir les résultats et écrire dans le fichier
for (i in seq_along(results)) {
  print(i)
  if (i == 6) {
    cat(paste(results[[i]]$combination, collapse = ","), "\t", file = output_file)
    cat(length(results[[i]]$intersection), "\n", file = output_file)
    break
  }
  for (k in seq_along(results[[i]])) {
    cat(paste(results[[i]][[k]]$combination, collapse = ","), "\t", file = output_file)
    cat(length(results[[i]][[k]]$intersection), "\n", file = output_file)
  }
  if (i < length(results)) {
    cat("\n", file = output_file)
  }
}

# Fermer le fichier
close(output_file)

data <- read.delim(gsub('XXX',output,'Shared_credible_SNPs_XXX.tsv'))

# Ajouter une colonne indiquant le nombre d'éléments comparés dans chaque combinaison
data$Num_Comparisons <- sapply(strsplit(as.character(data$Comparison), ","), length)


ggplot(data, aes(x = Comparison, y = Nbr_SNPs, fill = Comparison)) +
  geom_bar(stat = "identity") +
  labs(x = "Traits", y = "Number of SNPs", title = "Number of credible SNPs Shared between Different Trait") +
  theme_minimal()+
  scale_x_discrete(labels = function(x) gsub(",", "\n", x)) +
  facet_wrap(~ Num_Comparisons, ncol = 2, scales = "free")+
  theme(legend.position = "none", axis.text.x = element_text(size = 7))

ggsave(gsub('XXX',output,'Shared_credible_SNPs_XXX_plot.pdf'))



#### Ecrire fichier qui contient ID_SNP Nbr_Trait pour chaque SNPs credible dans au moins un dataframe
# Créer une liste pour stocker le compte de chaque SNP
snp_count <- list()

# Parcourir chaque fichier
for (i in 1:length(dfs)) {
  # Extraire les SNP crédibles du fichier
  snps <- dfs[[i]][dfs[[i]]$credible == "Credible", "X"]
  # Mettre à jour le compte de chaque SNP
  for (snp in snps) {
    if (is.null(snp_count[[snp]])) {
      snp_count[[snp]] <- 1
    } else {
      snp_count[[snp]] <- snp_count[[snp]] + 1
    }
  }
}

# Créer un dataframe avec les SNP et leur compte
snp_df <- data.frame(
  SNP = names(snp_count),
  Count = unlist(snp_count)
)
# Diviser la colonne SNP_CHROM en deux colonnes distinctes pour le SNP et le CHROM
snp_df <- separate(snp_df, SNP, into = c("CHROM","POS"), sep = "_", remove = FALSE)


# Écrire le dataframe dans un fichier CSV
write.csv(snp_df, gsub('XXX',output,"SNP_Count_XXX.csv"),quote=FALSE, row.names = FALSE)

all_snps <- dfs[[1]][,'X']
not_credible <- data.frame(SNP=all_snps[!(all_snps %in% snp_df$SNP)],Count=0)
not_credible  <- separate(not_credible , SNP, into = c("CHROM","POS"), sep = "_", remove = FALSE)

SNP_df <- data.frame(rbind(not_credible,snp_df))
# a. Manhattan plot bayésien

# Convertir la variable POS en numérique
SNP_df$POS <- as.numeric(SNP_df$POS)

# Créer le graphique en utilisant la palette de couleurs
palette <- rainbow(length(unique(SNP_df$Count)))

repartition <- ggplot(SNP_df, aes(x = POS, y = Count, color = factor(Count))) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Repartition SNPs along Chromosomes",
       x = "SNPs",
       y = "Number of trait associated") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        axis.ticks.x = element_blank()) +
  facet_wrap(~ CHROM, ncol = 2, scales = "free_x") +
  scale_x_continuous(breaks = seq(0, max(SNP_df$POS), by = 1000000))+
  theme(legend.position = "none")

print(repartition)

ggsave(gsub('XXX',output,'Repartition_SNPs_along_K_XXX.pdf'))

          