### Script to compare credible SNP between MT and ST
output <- 'VanRaden_A6_NGM_0.99'


library(ggplot2)
library(matrixStats)
library(tidyr)
library(gridExtra)
library(grid)


setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned")
#setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")
# Charger les fichiers CSV
files <- list.files(pattern = gsub('XXX',output,"Summary_[A-Z]*_XXX.csv"))
file_outputs_MT <- sub("Summary_(.*)_VanRaden_A6_NGM_0.99.csv", "\\1", files)
dfs_MT <- lapply(files, read.csv)



setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")
#setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
files <- list.files(pattern = gsub('XXX',output,"Summary_[A-Z]*_Uni_XXX.csv"))
file_outputs_ST <- sub("Summary_(.*)_Uni_VanRaden_A6_NGM_0.99.csv", "\\1", files)
dfs_ST <- lapply(files, read.csv)

setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
## Keep only SNP that are only associate with one trait
snp_df <- read.csv(gsub('XXX',output,"SNP_Count_Uni_XXX.csv"))
snp_df <- snp_df[snp_df$Count==1,]

setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")
snp_df_MT <- read.csv(gsub('XXX',output,"SNP_Count_XXX.csv"))

setwd("~/Documents/Worms/GBLUP")



# Combine the MT and ST data frames
combined_dfs <- list()
for (i in 1:length(dfs_MT)) {
  df_MT <- dfs_MT[[i]]
  df_ST <- dfs_ST[[i]]
  
  df_MT$source <- "MT"
  df_ST$source <- "ST"
  
  combined_df <- rbind(df_MT, df_ST)
  combined_dfs[[i]] <- combined_df
}


# Plotting
plots <- list()
plots_mp <- list()

for (i in 1:length(combined_dfs)) {
  df <- combined_dfs[[i]]
  df <- separate(df, X, into = c("CHROM", "POS"), sep = "_", remove = FALSE)
  df$POS <- as.numeric(df$POS)
  df$CHROM <- gsub('23','X',df$CHROM)
  
  trait <- file_outputs_MT[[i]]
  
  # Add new credible status for coloring
  df$credible_status <- "Not Credible"
  df$credible_status[df$credible == 'Credible' & df$source == 'MT'] <- "Credible MT"
  df$credible_status[df$credible == 'Credible' & df$source == 'ST'] <- "Credible ST"
  df$credible_status <- factor(df$credible_status, levels = c("Not Credible", "Credible MT", "Credible ST"))
  
  # Manhattan plot for combined data
  manhattan_plot <- ggplot(df, aes(x = POS, y = median, color = credible_status, alpha = credible_status)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = c("Not Credible" = 'grey', "Credible MT" = 'red', "Credible ST" = 'blue')) +
    scale_alpha_manual(values = c("Not Credible" = 0.2, "Credible MT" = 1, "Credible ST" = 1)) +
    theme_minimal() +
    ggtitle(paste("Manhattan Plot Bayésien", trait)) +
    labs(x = "SNPs",
         y = "Effet médian à posteriori") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    facet_wrap(~ CHROM, ncol = 2, scales = "free_x")
  
  ggsave(paste0('Manhattan_plots_', output, '-',trait,'.pdf'))
  plots_mp[[i]] <- manhattan_plot
  
  # Density plot for credible SNPs
  df_cred_MT <- df[df$credible_status == 'Credible MT',]
  df_cred_ST <- df[df$credible_status == 'Credible ST',]
  
  g <- ggplot() +
    geom_density(data = df[df$credible_status == 'Not Credible',], aes(x = median, fill = 'Not Credible'), color = 'grey', alpha = 0.3) +
    geom_density(data = df_cred_MT, aes(x = median, fill = 'Credible MT'), color = 'red', alpha = 0.3) +
    geom_density(data = df_cred_ST, aes(x = median, fill = 'Credible ST'), color = 'blue', alpha = 0.3) +
    scale_fill_manual(name = "Status", values = c("Not Credible" = 'grey', "Credible MT" = 'red', "Credible ST" = 'blue')) +
    labs(title = paste("Density median SNPs effects -", trait),
         x = "Médiane",
         y = "Densité") +
    theme_minimal()
  
  plots[[i]] <- g
}

# Function to extract the legend from a ggplot
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Get the legend from the first plot
legend_mp <- get_legend(plots_mp[[1]])

# Remove the legend from all plots
plots_mp <- lapply(plots_mp, function(plot) plot + theme(legend.position = "none"))

# Arrange plots with the legend
pdf(paste0('Manhattan_plots_', output, '.pdf'), width = 14, height = 10)
grid.arrange(grobs = plots_mp, ncol = 2, top = textGrob("Manhattan Plots", gp = gpar(fontsize = 20, fontface = "bold")))
dev.off()

pdf('Legend_MP.pdf')
grid.arrange(legend_mp)
dev.off()

# Arrange density plots
# Get the legend from the first plot
legend <- get_legend(plots[[1]])

# Remove the legend from all plots
plots <- lapply(plots, function(plot) plot + theme(legend.position = "none"))

pdf(paste0('Density_plots_', output, '.pdf'), width = 14, height = 10)
grid.arrange(grobs = plots, ncol = 2, top = textGrob("Density Plots", gp = gpar(fontsize = 20, fontface = "bold")))
dev.off()

pdf('Legend_density.pdf')
grid.arrange(legend)
dev.off()


# Créer un dataframe combiné pour tous les bar plots
all_bar_data <- data.frame()

for (i in 1:length(combined_dfs)) {
  df <- combined_dfs[[i]]
  df <- separate(df, X, into = c("CHROM", "POS"), sep = "_", remove = FALSE)
  df$POS <- as.numeric(df$POS)
  
  df$credible_status <- "Not Credible"
  df$credible_status[df$credible == 'Credible' & df$source == 'MT'] <- "Credible MT"
  df$credible_status[df$credible == 'Credible' & df$source == 'ST'] <- "Credible ST"
  df$credible_status <- factor(df$credible_status, levels = c("Not Credible", "Credible MT", "Credible ST"))
  
  trait <- file_outputs_MT[[i]]
  
  # Compare SNPs credible in both conditions
  df_cred_MT <- df[df$credible_status == 'Credible MT',]
  df_cred_ST <- df[df$credible_status == 'Credible ST',]
  
  snps_cred_MT <- df_cred_MT$X
  snps_cred_ST <- df_cred_ST$X
  common_snps <- intersect(snps_cred_MT, snps_cred_ST)
  unique_MT <- setdiff(snps_cred_MT, snps_cred_ST)
  unique_ST <- setdiff(snps_cred_ST, snps_cred_MT)
  
  # SNPs only in ST but not shared by other traits
  unique_ST_not_shared <- unique_ST[!(unique_ST %in% snp_df$SNP)]
  # SNPs only in ST but shared by other traits
  unique_ST_shared <- unique_ST[unique_ST %in% snp_df$SNP]
  
  # Préparer les données pour le bar plot
  bar_data <- data.frame(
    Category = c("Only MT", "Common", "Only ST (Shared)", "Only ST (Not Shared)"),
    Count = c(length(unique_MT), length(common_snps), length(unique_ST_shared), length(unique_ST_not_shared)),
    Trait = trait
  )
  
  all_bar_data <- rbind(all_bar_data, bar_data)
}

# Créer le bar plot combiné
combined_bar_plot <- ggplot(all_bar_data, aes(x = Category, y = Count, fill = Trait)) +
  geom_bar(stat = "identity", position = "stack", colour = "black") +
  labs(title = "Comparison of Credible SNPs",
       x = "Category",
       y = "Count") +
  theme_minimal()

# Sauvegarder le bar plot combiné
pdf(paste0('Combined_Bar_plots_', output, '.pdf'), width = 14, height = 10)
print(combined_bar_plot)
dev.off()


# Loop through each dataframe in the combined_dfs list
for (i in 1:length(combined_dfs)) {
  # Extract the current dataframe
  df <- combined_dfs[[i]]
  
  # Separate the 'X' column into 'CHROM' and 'POS'
  df <- separate(df, X, into = c("CHROM", "POS"), sep = "_", remove = FALSE)
  df$POS <- as.numeric(df$POS)
  
  # Assign credible status based on 'credible' and 'source' columns
  df$credible_status <- "Not Credible"
  df$credible_status[df$credible == 'Credible' & df$source == 'MT'] <- "Credible MT"
  df$credible_status[df$credible == 'Credible' & df$source == 'ST'] <- "Credible ST"
  df$credible_status <- factor(df$credible_status, levels = c("Not Credible", "Credible MT", "Credible ST"))
  
  # Get the trait name from file_outputs_MT
  trait <- file_outputs_MT[[i]]
  
  # Filter dataframes for credible SNPs in MT and ST
  df_cred_MT <- df[df$credible_status == 'Credible MT',]
  df_cred_ST <- df[df$credible_status == 'Credible ST',]
  
  # Extract SNPs that are credible in MT and ST
  snps_cred_MT <- df_cred_MT$X
  snps_cred_ST <- df_cred_ST$X
  
  # Find common and unique SNPs between MT and ST
  common_snps <- intersect(snps_cred_MT, snps_cred_ST)
  unique_MT <- setdiff(snps_cred_MT, snps_cred_ST)
  unique_ST <- setdiff(snps_cred_ST, snps_cred_MT)
  
  # SNPs only in ST and not shared by other traits
  unique_ST_not_shared <- unique_ST[!(unique_ST %in% snp_df$SNP)]
  
  # SNPs only in ST but shared by other traits
  unique_ST_shared <- unique_ST[unique_ST %in% snp_df$SNP]
  
  # Subset snp_df_MT for SNPs in unique_ST_not_shared
  snp_df_MT_sub <- snp_df_MT[snp_df_MT$SNP %in% unique_ST_not_shared,]
  
  # Print the results
  cat(paste(
    trait, 
    'Number of traits in ST_not_shared:', length(unique_ST_not_shared), 
    '\nNumber of them pleiotropic in MT (associated with another trait):', sum(snp_df_MT_sub$Count != 1), 
    '\nNon-pleiotropic in MT (associated with another trait):',sum(snp_df_MT_sub$Count == 1),
    # '\nNon-pleiotropic in MT (associated with another trait):', paste(snp_df_MT_sub[snp_df_MT_sub$Count == 1,]$SNP, collapse = ", "), 
    '\n\n'
  ))
}
