setwd("~/Documents/Worms/GBLUP")
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP")
matrix <- 'Inverted_kinship_matrix_VanRaden_A6.csv'
output <- 'A6_NaCl_NGM'

library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

Ginv <- as.matrix(read.csv(matrix))
colnames(Ginv) <- gsub('CeMee', '', colnames(Ginv))
colnames(Ginv) <- gsub('_sorted', '', colnames(Ginv))
rownames(Ginv) <- colnames(Ginv)

pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
pheno <- pheno[, c('pop_label', 'temperature', 'rel_humidity', "session_id",
                   'logD', 'date_str', 'env_label',
                   "SF", "SB", "FS", "FB", "BS", "BF")]

## Kept only lines with phenotypes and are in Ginv
pheno_subset <- pheno[pheno$pop_label %in% colnames(Ginv),]
pheno_subset <- pheno_subset[pheno_subset$env_label %in% c('NaCl', 'NGM'),]

# Scale covariates and reduce traits
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")

pheno_subset$temperature <- as.numeric(pheno_subset$temperature)
pheno_subset$rel_humidity <- as.numeric(pheno_subset$rel_humidity)
pheno_subset$logD <- as.numeric(pheno_subset$logD)

for (j in c('temperature', "rel_humidity", "logD")) {
  pheno_subset[, j][is.na(pheno_subset[, j])] <- mean(pheno_subset[, j], na.rm = TRUE)
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

for (j in vect_P_traits) {
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

# Plot phenotype data
pheno_long <- pheno_subset %>%
  pivot_longer(cols = c(FS, SF, BS, SB, BF, FB), names_to = "variable", values_to = "value")

pheno_long$variable <- factor(pheno_long$variable, levels = c("FS", "SF", "BS", "SB", "BF", "FB"))

# Créer le graphique avec ggplot2
g <- ggplot(pheno_long, aes(x = value, fill = env_label, color = env_label)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scales = "free", ncol = 3) +  # Utiliser ncol pour organiser les facettes
  theme_bw() +
  labs(title = "Trait Distribution by Condition",
       x = "Value",
       y = "Density") +
  theme(legend.position = "bottom")  # Placer la légende en bas

ggsave(gsub('XXX', output, 'XXX_plots_pheno_subset.pdf'))

# Calculer la moyenne
pheno_long_mean <- pheno_long %>%
  group_by(pop_label, variable) %>%
  summarize(value = mean(value, na.rm = TRUE))

# Initialiser le DataFrame final
df_to_plot <- data.frame(Line = unique(pheno_long_mean$pop_label),
                         SF = rep(0, length(unique(pheno_long_mean$pop_label))),
                         SB = rep(0, length(unique(pheno_long_mean$pop_label))),
                         FS = rep(0, length(unique(pheno_long_mean$pop_label))),
                         FB = rep(0, length(unique(pheno_long_mean$pop_label))),
                         BS = rep(0, length(unique(pheno_long_mean$pop_label))),
                         BF = rep(0, length(unique(pheno_long_mean$pop_label))))

# Remplir le DataFrame avec les valeurs moyennes
for(i in vect_P_traits){
  df_to_plot[,i] <- pheno_long_mean[pheno_long_mean$variable == i,]$value
}
comb<-combn(vect_P_traits,2)
plots <- list()

# Using the Set3 color palette from RColorBrewer
base_colors <- brewer.pal(9, "Set1")

#Augmenting the palette with additional colors
colors <- c(base_colors, rainbow(15 - length(base_colors)))

for(i in 1:ncol(comb)){
  x <- comb[, i][[1]]
  y <- comb[, i][[2]]
  
  corr_coef <- cor(df_to_plot[[x]], df_to_plot[[y]], method = "pearson") # ou method = "spearman"
  
  p <- ggplot(df_to_plot, aes_string(x = x, y = y, color = factor(i))) +
    geom_point(show.legend = FALSE) +
    labs(x = x, y = y) +
    scale_color_manual(values = colors[i]) + 
    ggtitle(paste(x, "vs", y),subtitle = paste("Correlation:", round(corr_coef, 2))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  plots[[i]] <- p
  
  
  plots[[i]] <- p
}

pdf(gsub('XXX',output,"Pheno_Trait_VS_traits_plot_XXX.pdf"),width = 9)
grid.arrange(grobs = plots, ncol = 5)
dev.off()


  
