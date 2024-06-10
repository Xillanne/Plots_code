## Script which compare credible lines for BV

file<-'Lines_credibles_VanRaden_A6_NGM.csv'

library(dplyr)

setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
ST <- read.csv(file)
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")
MT <- read.csv(file)

cat(paste('Number of line credibled for a specific trait shared in MT and ST :',nrow(semi_join(ST,MT,by=c('Line','Trait')))))
cat(paste('Number of line credibled for a specific trait MT and not in ST :',nrow(anti_join(MT,ST,by=c('Line','Trait')))))
cat(paste('Number of line credibled for a specific trait ST and not in MT :',nrow(anti_join(ST,MT,by=c('Line','Trait')))))