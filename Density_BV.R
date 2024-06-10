#### RESULTS #### https://people.bath.ac.uk/jjf23/mixchange/
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(HDInterval)
library(combinat)
library(RColorBrewer)

output <- 'VanRaden_A6_NGM'
Sol_file <- 'VanRaden_A6_NGM_MCMCmodel_Sol.csv'
M_file <- 'A6_t_convert_genotype.csv'

Sol <- list()
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")
Sol[[1]] <- read.csv(Sol_file)
Sol_name <- c('Multi')
M <- as.data.frame(fread(M_file))
rownames(M) <- M[,1]
M[,1] <- NULL
M <- as.matrix(M)

setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
vect_P_traits <- c("SF",
                   "SB",
                   "FS",
                   "FB",
                   "BS",
                   "BF")
for(p in vect_P_traits){
  Sol_name <- c(Sol_name,'Uni')
  Sol[[length(Sol)+1]] <- read.csv(paste0(output,'_',p,'_MCMCmodel_Sol.csv'))
}

pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
M <- M[rownames(M) %in% pheno$pop_label,] #keep only the name of line which were phenotyped and genotyped

setwd("~/Documents/Worms/GBLUP")
### Trait breeding value mean for each line ###
hist_breeding_values <- data.frame(Line=character(0),
                                   Trait=character(0),
                                   Median=numeric(0),
                                   Var=numeric(0),
                                   Lower_ic = numeric(0),
                                   Upper_ic = numeric(0),
                                   Credible = character(0),
                                   source=character(0)
                                   )
is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}


i <- 1
for(model_MCMC_WI_Sol in Sol){
  for( name in rownames(M)){
    ## Print distribution for all Lines 
    # Select columns which contains name of the line (ie breeding value distribution for each trait)
    col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(name, "$"),colnames(model_MCMC_WI_Sol))] 
    if(length(col_subset)!=0){ # if the line is contains
      subset <- data.frame(model_MCMC_WI_Sol[,col_subset]) #consider just this line
      colnames(subset)<-gsub(paste0("^", "trait", "|", ".pop_label.",name, "$"), "", colnames(subset)) #rename with FS SF BF..
      
      
      ## Find line with a credible breeding value
      # Calculate credible interval
      HDI <- hdi(subset) #95%
      
      if(length(sub('trait([A-Z]*).pop_label.*','\\1',col_subset))!=1){
        trait <- sub('trait([A-Z]*).pop_label.*','\\1',col_subset)
      } else {
        trait <- vect_P_traits[[i-1]]
      }
      # Compute a summary for the line
      summary_df <- data.frame(
        Line = rep(name,length(col_subset)),
        Trait = trait,
        Median = colMedians(as.matrix(subset)),
        Sd = colSds(as.matrix(subset)),
        Lower_ic = HDI["lower",],
        Upper_ic = HDI["upper",],
        Credible = !apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"])), #TRUE if 0 not in interval
        source=rep(Sol_name[i],length(col_subset))
      )
      
      # Add it with other
      hist_breeding_values <- data.frame(rbind(hist_breeding_values,summary_df))
      }
  }
  i <- i +1
}

# Définir les couleurs pour chaque catégorie avec les nouveaux labels
colors <- c("MT All" = "blue", "MT Credible" = "cyan", "ST All" = "pink", "ST Credible" = "red")

# Créer le graphique avec normalisation des densités et les nouveaux labels
m<-ggplot() +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Multi',], aes(x = Median, fill = 'MT All'), color = 'black', alpha = 0.4) +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Uni',], aes(x = Median, fill = 'ST All'), color = 'black', alpha = 0.4) +
  labs(title = "Density of Median Breeding values",
       x = "Median",
       y = "Density",
       fill = "Legend") +
  scale_fill_manual(values = colors) +
  theme_minimal()


# Créer le graphique avec normalisation des densités et les nouveaux labels
s<-ggplot() +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Multi',], aes(x = Sd, fill = 'MT All'), color = 'black', alpha = 0.4) +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Uni',], aes(x = Sd, fill = 'ST All'), color = 'black', alpha = 0.4) +
  labs(title = "Density of Sd breeding values",
       x = "Sd",
       y = "Density",
       fill = "Legend") +
  scale_fill_manual(values = colors) +
  theme_minimal()

ggsave(gsub('XXX',output,'DensityPlot_BV_XXX.pdf'))

grid.arrange(m,s,nrow=1)

