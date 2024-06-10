#### RESULTS #### https://people.bath.ac.uk/jjf23/mixchange/
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(HDInterval)
library(combinat)
library(RColorBrewer)

setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")

Sol_file <- 'VanRaden_A6_NGM_MCMCmodel_Sol.csv'
M_file <- 'A6_t_convert_genotype.csv'
output <- 'VanRaden_A6_NGM'

model_MCMC_WI_Sol <- read.csv(Sol_file)
M <- as.data.frame(fread(M_file))
rownames(M) <- M[,1]
M[,1] <- NULL
M <- as.matrix(M)

pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
M <- M[rownames(M) %in% pheno$pop_label,] #keep only the name of line which were phenotyped and genotyped

### Trait breeding value mean for each line ###
hist_breeding_values <- data.frame(Line=character(0),
                                   Trait=character(0),
                                   Median=numeric(0),
                                   SD=numeric(0))
distrib <- list() #Will contains distribution of breeding value for each trait per line

is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}

vect_P_traits <- c("SF",
                   "SB",
                   "FS",
                   "FB",
                   "BS",
                   "BF")

i<-1
for( name in rownames(M)){
  ## Print distribution for all Lines 
  # Select columns which contains name of the line (ie breeding value distribution for each trait)
  col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(name, "$"),colnames(model_MCMC_WI_Sol))] 
  if(length(col_subset)!=0){ # if the line is contains
    subset <- data.frame(model_MCMC_WI_Sol[,col_subset]) #consider just this line
    colnames(subset)<-gsub(paste0("^", "trait", "|", ".pop_label.",name, "$"), "", colnames(subset)) #rename with FS SF BF..
    
    
    sub_melt <- melt(subset) #Transform in two columns dataframe
    colnames(sub_melt)<-c('traits','breeding_value') #rename columns
    # Plot posterior distribution for each trait for specific line
    g<-ggplot(sub_melt,aes(x = breeding_value,color=traits))+
      geom_density()+
      ggtitle(paste0('Distribution of breeding value for each trait for ',name)) + 
      theme_bw()
    # Save the plot with name of the line
    distrib[[i]] <- list(line=name,
                         plot=g)
    i<-i+1
    
    
    ## Find line with a credible breeding value
    # Calculate credible interval
    HDI <- hdi(subset) #95%
    
    # Compute a summary for the line
    summary_df <- data.frame(
      Line = rep(name,6),
      Trait = vect_P_traits,
      Median = colMedians(as.matrix(subset)),
      SD = colSds(as.matrix(subset)),
      Lower_ic = HDI["lower",],
      Upper_ic = HDI["upper",],
      Credible = !apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"])) #TRUE if 0 not in interval
    )
    
    # Add it with other
    hist_breeding_values <- data.frame(rbind(hist_breeding_values,summary_df))
    }
}


# Plot hist_breeding_values only for credible lines
Median_trait <- ggplot(hist_breeding_values, aes(x = Line, y = Median, color = Credible, group = Trait)) +
  geom_point(position = position_dodge(width = 0.9), size = 2) +
  geom_errorbar(aes(ymin = Median - SD, ymax = Median + SD), width = 0.2, linewidth=0.25, position = position_dodge(width = 0.9), color = "black") +
  labs(x = "Lines", y = "Median", color = "Credible") +
  theme_minimal() +
  theme(axis.text.x = element_blank())+  # Remove x-axis text) +
  facet_wrap(~Trait, ncol = 2, nrow = 3)+
  ggtitle(gsub('XXX',output,'Breeding values XXX'))

Median_trait
pdf_filename <- gsub('XXX',output,"Trait_distribution_each_lines_XXX.pdf")
ggsave(pdf_filename, Median_trait)


##### Forest plot for credible line
cred_BV <- hist_breeding_values[hist_breeding_values$Credible==TRUE,]
write.csv(cred_BV,gsub('XXX',output,'Lines_credibles_XXX.csv'),row.names = FALSE,quote = FALSE)
length(unique(cred_BV$Line))

# Définir les couleurs pour chaque Trait
colors <- c("SF" = "blue", "FS" = "green", "FB" = "red","SB" = "cyan", "BS" = "pink", "BF" = "orange")

ggplot(cred_BV, aes(x = Median, y = Line, color = Trait)) +
  geom_point(position = position_dodge(width = 0.25), size = 2) +  # Points pour les moyennes avec décalage
  #geom_errorbarh(aes(xmin = Lower_ic, xmax = Upper_ic), size=0.5,height = 0.2, position = position_dodge(width = 0.25)) +  # Barres d'erreur horizontales avec décalage
  #geom_segment(aes(x = -0.7, xend = 0.7, y = Line, yend = Line), color = "black", alpha = 0.5) +  # Tracer une ligne entre chaque y
  scale_color_manual(values = colors) +  # Appliquer les couleurs définies
  labs(title = "Forest Plot des Traits par Ligne",
       x = "Valeur Moyenne",
       y = "Ligne",
       color = "Trait") +
  theme_minimal()+  # Utiliser un thème minimal
  theme(panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_line(color = "gray"))
  
ggsave(gsub('XXX',output,"forest_plot_XXX.pdf"), width = 4, height = 8)


### Plot to compare breeding value per trait : for each line plot the value in a trait against an other
df_to_plot <- data.frame(Line=unique(hist_breeding_values$Line),
                         SF=rep(0,nrow(hist_breeding_values)),
                         SB=rep(0,nrow(hist_breeding_values)),
                         FS=rep(0,nrow(hist_breeding_values)),
                         FB=rep(0,nrow(hist_breeding_values)),
                         BS=rep(0,nrow(hist_breeding_values)),
                         BF=rep(0,nrow(hist_breeding_values)))

to_compare <- unique(hist_breeding_values$Trait)
for(i in to_compare){
  df_to_plot[,i] <- hist_breeding_values[hist_breeding_values$Trait==i,]$Median
}
comb<-combn(to_compare,2)
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
    ggtitle(paste(x, "vs", y), subtitle = (paste("Correlation:", round(corr_coef, 2)))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  plots[[i]] <- p
  
  
  plots[[i]] <- p
}

pdf(gsub('XXX',output,"Trait_VS_traits_plot_XXX.pdf"),width = 9)
grid.arrange(grobs = plots, ncol = 5)
dev.off()
