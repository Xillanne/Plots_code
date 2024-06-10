setwd("~/Documents/Worms/GBLUP/G-matrix")
name <- 'A6_NaCl'

to_seek <- gsub('XXX',name,".*_XXX_G_matrix.csv")
files <- list.files(pattern = to_seek)
files_name <- sub(paste0('(.*)_',name,"_G_matrix.csv"), "\\1", files)
G_kin_matrices <- lapply(files, read.csv)

Pheno <- read.csv(gsub('XXX',name,"XXX_MCMCmodel_VCV_G_matrix.csv"))


library(ggplot2)
library(reshape2)
library(combinat)
library(gridExtra)


plots <- list() 

for(i in 1:length(files)){
  df1_melt <- melt(G_kin_matrices[[i]], varnames = c("Row", "Col"), value.name = "Value1")
  
  df2_melt <- melt(Pheno, varnames = c("Row", "Col"), value.name = "Value2")
  
  # Combine the two melted data frames by row and column indices
  combined_df <- data.frame(cbind(df1_melt,df2_melt))
  
  # Plot the values of matrix1 against matrix2
  g <- ggplot(combined_df, aes(x = Value1, y = Value2)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Add diagonal line
    ggtitle(paste0(comb_name[,i][1],' VS Phenotypic')) +
    xlab(files_name[i]) +
    ylab('Phenotypic') +
    theme_bw()
  
  plots[[i]] <- g
}

pdf(gsub('XXX',name,"Comparison_G_matrix_XXX.pdf"))
grid.arrange(grobs = plots, ncol = 2)
dev.off()
