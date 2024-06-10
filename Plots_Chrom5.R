## Script to plot for all pop zoom in chrom 5
setwd("~/Documents/Worms/VariantCalling/Hetero")

library(tidyr)
library(gridExtra)

# Charger les fichiers CSV
files <- list.files(pattern = "Heterozygosity_.*.csv")
populations <- sub("Heterozygosity_(.*).csv",'\\1',files)
dfs <- lapply(files, read.csv)

## Load HD regions
df <- read.delim("20220216_c_elegans_divergent_regions_strain.bed", header = FALSE)
df$V4 <- NULL
df <- unique(df)

# Function to merge overlapping intervals
merge_intervals <- function(df) {
  result <- list()
  current_interval <- c(df$V2[1], df$V3[1])
  
  for (i in 2:nrow(df)) {
    if (df$V2[i] <= current_interval[2]) {
      current_interval[2] <- max(current_interval[2], df$V3[i])
    } else {
      result <- rbind(result, current_interval)
      current_interval <- c(df$V2[i], df$V3[i])
    }
  }
  
  result <- rbind(result, current_interval)
  return(as.data.frame(result))
}

# Store chrom5 results for each populations
Df_chrom5 <- data.frame(X=integer(0),
                        Chromosome = character(0),
                        Bin_number = integer(0),
                        Start = integer(0),
                        End = integer(0),
                        Number_of_heterozygotes = numeric(0),
                        HD = numeric(0),
                        Bin_chrom = integer(0),
                        Population = character(0))

plots <- list()

for(f in 1:length(files)){
  pop <- populations[[f]]
  Results <- dfs[[f]]
  Results[is.na(Results)] <- 0
  
  # Apply function to merge intervals
  chrom <- unique(Results$Chromosome)
  
  HD_df <- data.frame(Chromosome = character(0), Start = integer(0), End = integer(0))
  for (i in chrom) {
    df_chrom <- df[df$V1 == i, ]
    res <- merge_intervals(df_chrom)
    colnames(res) <- c("Start", "End")
    res$Chromosome <- rep(i, nrow(res))
    HD_df <- rbind(HD_df, res)
  }
  
  row.names(HD_df) <- NULL
  
  # Initialize HD_r vector with the length of Results$Chromosome
  HD_r <- integer(length(Results$Chromosome))
  Bin_chrom <- integer(length(Results$Chromosome))
  # Initialize nbr variable to keep track of the number of rows processed
  nbr <- 0
  
  # Loop over each chromosome in the chrom vector
  for (c in chrom) {
    # Subset Results and HD_df data frames based on the current chromosome
    Results_chrom <- Results[Results$Chromosome == c, ]
    HD_df_chrom <- HD_df[HD_df$Chromosome == c, ]
    
    # Loop over each row of HD_df_chrom to iterate through the regions
    for (j in 1:nrow(HD_df_chrom)) {
      # Extract start and end coordinates of the current region
      start <- HD_df_chrom$Start[[j]]
      end <- HD_df_chrom$End[[j]]
      
      # Loop over each row of Results_chrom to check for overlap with the current region
      for (i in 1:nrow(Results_chrom)) {
        Bin_chrom[i+nbr] <- i
        #print(c(start,end))
        #print(c(Results_chrom$End[i],Results_chrom$Start[i]))
        
        # Check if there is overlap between the current region and the current result
        if (start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
            start <= Results_chrom$Start[i] && end >= Results_chrom$End[i] || 
            start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
            start >= Results_chrom$Start[i] && end <= Results_chrom$End[i]) {
          # If there is overlap, mark the corresponding entry in HD_r as -1
          #print(i + nbr)
          HD_r[i + nbr] <- -0.05
          res <- c(res,1)
        }
        #if(Results_chrom$End[i]==1007846){break}
      }
      #print(i)
      #break
    }
    #print(head(HD_r,n = 300))
    # Update nbr to keep track of the total number of rows processed
    nbr <- nbr + nrow(Results_chrom)
  }
  HD_r[HD_r==0] <- NaN #Set at missing in odrer to not plot them
  Results$HD <- HD_r
  Results$Bin_chrom <- Bin_chrom
  
  # Convert Chromosome column to a factor if it isn't already
  Results$Chromosome <- as.factor(Results$Chromosome)
  
  Chrom5 <- Results[Results$Chromosome=='V',]
  Chrom5_zoom <- Chrom5
  #Chrom5_zoom <- Chrom5[Chrom5$Start>16813367,]
  #Chrom5_zoom <- Chrom5_zoom[Chrom5_zoom$End<17627977,]
  
  Chrom5$Population <- rep(pop,nrow(Chrom5))
  Df_chrom5 <- data.frame(rbind(Df_chrom5,Chrom5))
  
  g5<-ggplot(Chrom5_zoom)+
    geom_point(aes(x=Start,y=Number_of_heterozygotes), alpha = 0.6,color='blue')+
    geom_point(aes(x = Start, y = HD), color = "black", alpha = 0.6) +
    ggtitle(pop) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  plots[[f]] <- g5
}


ggplot(Df_chrom5)+
  geom_point(aes(x=Start,y=Number_of_heterozygotes,color=Population), alpha = 0.6) +
  geom_point(aes(x = Start, y = HD), color = "black", alpha = 0.6)+
  facet_wrap(~ Population, ncol = 3)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

ggsave('Chromosome_5_Allpop.pdf')
