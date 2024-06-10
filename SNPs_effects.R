### Script to analyse SNP effects
output <- 'VanRaden_A6_NaCl_0.99'
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned")
#setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")

library(ggplot2)
library(matrixStats)
library(tidyr)
library(gridExtra)



# Charger les fichiers CSV
files <- list.files(pattern = gsub('XXX',output,"Summary_[A-Z]*_XXX.csv"))
file_outputs <- sub("Summary_(.*)_VanRaden_A6_NaCl_0.99.csv", "\\1", files)
dfs <- lapply(files, read.csv)
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")

plots <- list()

## Keep only SNP that are only associate with this trait
snp_df <- read.csv(gsub('XXX',output,"SNP_Count_XXX.csv"))
snp_df <- snp_df[snp_df$Count==1,]

for(i in 1:length(dfs)){
  df <- dfs[[i]]
  df$credible <-factor(df$credible, levels = c("Not Credible", "Credible"))
  df <- separate(df , X, into = c("CHROM","POS"), sep = "_", remove = FALSE)
  df$POS <- as.numeric(df$POS)
  
  trait <- file_outputs[[i]]
  
  
  df_1 <- df[df$X %in% snp_df$SNP,]
  # a. Manhattan plot bayésien
  manhattan_plot <- ggplot(df_1, aes(x = POS, y = median, color = credible, alpha = credible)) +
    geom_point(size=0.5) + 
    # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    scale_color_manual(values = c("Not Credible" = 'grey', "Credible" = "red")) +
    scale_alpha_manual(values = c("Not Credible" = 0.01, "Credible" = 1)) +
    theme_minimal() +
    ggtitle(gsub('TTT',trait,"Manhattan Plot Bayésien TTT"))+
    labs(x = "SNPs",
         y = "Effet médian à posteriori") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    facet_wrap(~ CHROM, ncol = 2, scales = "free_x")
  #print(manhattan_plot)
  
  ggsave(gsub('TTT',trait,gsub('XXX',output,'Manhattan_plot_XXX_TTT.pdf')))
  
  df_cred <- df[df$credible=='Credible',]
  
  g <- ggplot() +
    geom_density(data = df, aes(x = median), color = 'blue', fill = 'blue', alpha = 0.3) +
    geom_density(data = df_cred, aes(x = median), color = 'red', fill = 'red', alpha = 0.3) +
    labs(title = paste("Density median SNPs effects -",trait),
         x = "Médiane",
         y = "Densité") +
    theme_minimal()
  plots[[i]] <- g
  
  
}

pdf(gsub('XXX',output,'Density_plot_XXX.pdf'))
grid.arrange(grobs=plots)
dev.off()

