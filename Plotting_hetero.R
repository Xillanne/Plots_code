args <- commandArgs(trailingOnly = TRUE)
Results_file <-  args[1] 
pop <- args[2]

setwd("~/Documents/Worms/VariantCalling/Hetero")
Results_file_1 <-  'Heterozygosity_GA.csv'
Results <- read.csv(Results_file_1)
Results_file <-  'Results_I_II_III_X/Heterozygosity_A6.csv'
pop <- 'A6'

Results_1 <- read.csv(Results_file)

Results <- data.frame(rbind(Results,Results_1))
Results$Bin_Number <- c(1:nrow(Results))
library(ggplot2)



Results[is.na(Results)] <- 0

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
        HD_r[i + nbr] <- -0.005
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

# Create the plot
g <- ggplot(data = Results) +
  geom_point(aes(x = Bin_chrom, y = Number_of_heterozygotes, color = Chromosome), alpha = 0.6) +
  geom_point(aes(x = Bin_chrom, y = HD), color = "black", alpha = 0.6) +
  ggtitle(paste0("Heterozygosity along chromosome for ", pop)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~ Chromosome, ncol = 2, scales = "free_x") # Adjust ncol to control the number of columns in the grid

# Print the plot
print(g)
pdf_filename <- paste0("Heterozygotity_HDregion_", gsub(" ", "_", pop), ".pdf")
ggsave(pdf_filename, g)


## Zoom on intresting regions
Chrom5 <- Results[Results$Chromosome=='V',]
Chrom5_zoom <- Chrom5[Chrom5$Start>16813367,]
Chrom5_zoom <- Chrom5_zoom[Chrom5_zoom$End<17627977,]
g5<-ggplot(Chrom5_zoom)+
  geom_point(aes(x=Start,y=Number_of_heterozygotes), alpha = 0.6,color='blue')+
  geom_point(aes(x = Start, y = HD), color = "black", alpha = 0.6) +
  ggtitle(paste0("Heterozygosity along chromosome for ", pop)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pdf_filename <- paste0("Heterozygotity_HDregion_Zoom_V_", gsub(" ", "_", pop), ".pdf")
ggsave(pdf_filename,g5)

library(ggplot2)

# Assuming 'Results' is your dataframe and it contains a column 'Chromosome'

# Replace NA values with 0
Results[is.na(Results)] <- 0

# Modify specific HD values as in the original code
Results[,Results$HD == -0.05] <- 1

# Define the number of repetitions for randomization
num_repetitions <- 50000

# Function for randomization
randomisation <- function(df, nbr) {
  ## Shuffle the heterozygosity and see what is the chance to be linked by chance
  df_random <- as.data.frame(df$HD)
  colnames(df_random) <- "HD"
  df_random$Number_of_heterozygotes <- sample(df$Number_of_heterozygotes) 
  df_sorted <- df_random[order(-df_random$Number_of_heterozygotes), ]
  rownames(df_sorted) <- c(1:length(df_sorted$Number_of_heterozygotes))
  df_sorted <- df_sorted[1:nbr, ]
  prop_hetero <- sum(-df_sorted$HD) / nbr
  return(prop_hetero)
}

# Initialize a list to store histograms
histograms <- list()

# Process each chromosome separately
for (chromosome in c('I','II','III','IV','V','X')) {
  # Filter the dataframe for the current chromosome
  chromosome_data <- Results[Results$Chromosome == chromosome, ]
  
  # Sort dataframe
  Results_sorted <- chromosome_data[order(-chromosome_data$Number_of_heterozygotes), ]
  rownames(Results_sorted) <- c(1:length(Results_sorted$Chromosome))
  
  # Get the 5% most important values
  nbr <- 5 * length(Results_sorted$Chromosome) / 100
  Results_sorted <- Results_sorted[1:nbr, ]
  prop_hetero <- sum(-Results_sorted$HD) / nbr
  
  # Apply the randomisation function and store the results in a vector
  distribution <- replicate(num_repetitions, randomisation(chromosome_data, nbr))
  
  # Create a histogram for the current chromosome
  hist <- ggplot(NULL, aes(x = distribution)) +
    geom_density(color = "black", fill = "white",bw=0.00015) +
    geom_vline(xintercept = prop_hetero, color = "red") +
    ggtitle(paste("Chromosome", chromosome)) +
    theme_bw()
  
  # Store the histogram in the list
  histograms[[chromosome]] <- hist
  
}


# Save the histogram to a PDF file
pdf(paste0("Randomisation_", pop, ".pdf"))
grid.arrange(grobs=histograms)
dev.off()


## Idem with the zoom
Chrom5_zoom[is.na(Chrom5_zoom)] <- 0
Chrom5_zoom[,Chrom5_zoom$HD==-0.05] <- 1

# Sort dataframe
Chrom5_zoom_sorted <- Chrom5_zoom[order(-Chrom5_zoom$Number_of_heterozygotes), ]
rownames(Chrom5_zoom_sorted) <- c(1:length(Chrom5_zoom$Chromosome))
# Get the 5% most important values
nbr <- 5 * length(Chrom5_zoom$Chromosome) / 100
Chrom5_zoom_sorted <- Chrom5_zoom_sorted[1:nbr, ]
#Chrom5_zoom_sorted[is.na(Chrom5_zoom_sorted)] <- 0
prop_hetero <- sum(-Chrom5_zoom_sorted$HD) / nbr


# Set the number of repetitions
num_repetitions <- 50000

# Apply the randomisation function 1000 times and store the Chrom5_zoom in a vector
distribution_Chrom5 <- replicate(num_repetitions, randomisation(Chrom5_zoom, nbr))

hist_Chrom5 <- ggplot(NULL, aes(x = distribution)) +
  geom_density(color = "black", fill = "white") +
  geom_vline(xintercept = prop_hetero, color = "red") +
  theme_bw()

print(hist_Chrom5 )
pdf_filename <- paste0("Randomisation_Chrom5_", gsub(" ", "_", pop), ".pdf")
ggsave(pdf_filename, hist_Chrom5 )


