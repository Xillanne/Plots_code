### Script to analyze results from Univariate GBLUP
library(MCMCglmm)

setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")

VCV_file <- 'VanRaden_A6_NaCl_BF_MCMCmodel_VCV.csv'

VCV <- read.csv(VCV_file)

posterior.mode(VCV)
