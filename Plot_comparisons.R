## Script to plot UNI VS MULTI

output <- 'VanRaden_A6_NaCl_0.99'

setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned")
#setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")
# Charger les fichiers CSV
files <- list.files(pattern = gsub('XXX',output,"Summary_[A-Z]*_XXX.csv"))
file_outputs <- sub("Summary_(.*)_VanRaden_A6_NaCl_0.99.csv", "\\1", files)
dfs_MT <- lapply(files, read.csv)



setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")
#setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
files <- list.files(pattern = gsub('XXX',output,"Summary_[A-Z]*_Uni_XXX.csv"))
file_outputs <- sub("Summary_(.*)_Uni_VanRaden_A6_NaCl_0.99.csv", "\\1", files)
dfs_ST <- lapply(files, read.csv)


setwd("~/Documents/Worms/GBLUP")

library(ggplot2)

# Combine the MT and ST data frames
combined_dfs <- data.frame()
combined_bis <- data.frame()
for (i in 1:length(dfs_MT)) {
  df_MT <- dfs_MT[[i]]
  df_ST <- dfs_ST[[i]]
  df_MT$IC <- abs(df_MT$lower - df_MT$upper)
  df_ST$IC <- abs(df_ST$lower - df_ST$upper)
  
  df_MT$source <- "MT"
  df_ST$source <- "ST"
  df_MT$Trait <- file_outputs[[i]]
  df_ST$Trait <- file_outputs[[i]]
  
  combined_df <- rbind(df_MT, df_ST)
  combined_dfs<- rbind(combined_df,combined_dfs)
  
  combined_bis <- rbind(combined_bis,merge(df_MT,df_ST,by = c('X','Trait','CHROM')))
}

names(combined_bis)[names(combined_bis) == "median.x"] <- "median_MT"
names(combined_bis)[names(combined_bis) == "median.y"] <- "median_ST"


# Définir les couleurs
colors <- c("ST All" = "blue", "ST Credible" = "cyan", "MT All" = "pink", "MT Credible" = "red")

# Créer le graphique avec les densités et les facettes
ggplot() +
  geom_density(data = combined_dfs[combined_dfs$source == 'MT',], aes(x = median, fill = 'MT All' ), color = 'black', alpha = 0.4) +
  geom_density(data = combined_dfs[combined_dfs$source == 'MT' & combined_dfs$credible == 'Credible',], aes(x = median, fill = 'MT Credible' ), color = 'black', alpha = 0.4) +
  geom_density(data = combined_dfs[combined_dfs$source == 'ST',], aes(x = median, fill = 'ST All' ), color = 'black', alpha = 0.4) +
  geom_density(data = combined_dfs[combined_dfs$source == 'ST' & combined_dfs$credible == 'Credible',], aes(x = median, fill = 'ST Credible' ), color = 'black', alpha = 0.4) +
  labs(title = "Density of Median SNP Effects",
       x = "Median",
       y = "Density",
       fill = "Legend") +
  facet_wrap(~Trait, ncol = 1) +
  scale_fill_manual(values = colors) +
  expand_limits(x = range(combined_dfs$median)) +
  theme_minimal()

ggsave(gsub('XXX',output,'Density_plots_comparions_Mt_Uni_XXX.pdf'))


# Ajouter une nouvelle colonne pour les couleurs
combined_bis$color <- with(combined_bis, ifelse(credible.y == 'Credible' & credible.x == 'Credible', 'Violet',
                                                ifelse(credible.y == 'Credible', 'Blue',
                                                       ifelse(credible.x == 'Credible', 'Red', 'Grey'))))

# Créer le graphique avec les points colorés en fonction de la nouvelle colonne
g <- ggplot(combined_bis, aes(x = median_ST, y = median_MT, color = color)) +
  geom_point(alpha = 0.6) +
  labs(title = "Plot of Median SNP Effects",
       x = "Median ST",
       y = "Median MT",
       color = "Credibility") +  # Titre de la légende
  facet_wrap(~Trait, ncol = 2) +  # Utiliser facet_wrap pour définir le nombre de colonnes
  theme_minimal() +
  scale_color_manual(values = c('Blue' = 'blue', 'Red' = 'red', 'Violet' = 'purple', 'Grey' = 'grey'),
                     labels = c('Blue' = 'Credible in ST', 'Red' = 'Credible in MT', 'Violet' = 'Credible in Both', 'Grey' = 'Not Credible')) +
  theme(legend.position = 'right')+  # Afficher la légende à droite
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")



ggsave(gsub('XXX',output,'Uni_VS_Multi_XXX.pdf'))

