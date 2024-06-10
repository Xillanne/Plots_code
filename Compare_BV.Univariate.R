## SINGLE TRAIT
#### RESULTS #### https://people.bath.ac.uk/jjf23/mixchange/
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(HDInterval)
library(combinat)

setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")

name <- 'VanRaden_A6_(0.[0-9]*)_NGM'
name1 <- 'VanRaden_A6_NGM'
M_file <- 'A6_t_convert_genotype.csv'

to_seek <- gsub('XXX',name,"XXX_MCMCmodel_Sol.csv")
files <- list.files(pattern = to_seek)
files_name <- sub(to_seek, "\\1", files)
dfs <- lapply(files, read.csv)
dfs[[length(files)+1]] <- as.data.frame(fread(gsub('XXX',name1,"XXX_MCMCmodel_Sol.csv")))
files_name<-c(files_name,'N')

M <- as.data.frame(fread(M_file))
rownames(M) <- M[,1]
M[,1] <- NULL
M <- as.matrix(M)

pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
M <- M[rownames(M) %in% pheno$pop_label,] #keep only the name of line which were phenotyped and genotyped

### Trait breeding value mean for each line ###
hist_breeding_values <- data.frame(Line=character(0),
                                   Pruned=character(0),
                                   Trait=character(0),
                                   Median=numeric(0),
                                   SD=numeric(0),
                                   Lower_ic = numeric(0),
                                   Upper_ic = numeric(0),
                                   Credible = character(0))
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

for(k in 1:length(files_name)){
  model_MCMC_WI_Sol <- dfs[[k]]
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
      i<-i+1
      
      
      ## Find line with a credible breeding value
      # Calculate credible interval
      HDI <- hdi(subset) #95%
      
      # Compute a summary for the line
      summary_df <- data.frame(
        Line = rep(name,6),
        Pruned = files_name[[k]],
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
}

ggplot(hist_breeding_values, aes(x=Pruned,y = Median, fill = Trait)) +
  geom_bar(stat = "identity") +
  labs(title = "Trait median along pruning",
       x = "Pruning",
       y = "Median") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Trait, ncol = 2)+
  theme(legend.position = "none")

ggsave('Median_BV_among_pruning.pdf')
