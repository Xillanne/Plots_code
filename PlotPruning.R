## Script to plot pruning results

setwd("~/Documents/Worms/GBLUP/Results_Pruning")

table <- read.csv2('summary_imputation.csv')

library(ggplot2)
library(reshape2)
colnames(table) <- gsub('nSNPS_afterF', 'P_1',colnames(table))

table_melt <- melt(table,value.name = 'Population')
colnames(table_melt) <- c('Population','Pruning','SNPs')

desired_order <- c('P_1', 'P_0.99', 'P_0.9','P_0.8','P_0.7','P_0.6','P_0.5','P_0.4','P_0.3','P_0.2','P_0.1') # replace with actual levels in desired order
table_melt$Pruning <- factor(table_melt$Pruning, levels = desired_order)

desired_order <- c('Founders','A6','CA','GA','GT','GM','LR','SMR') # replace with actual levels in desired order
table_melt$Population <- factor(table_melt$Population, levels = desired_order)

ggplot(table_melt, aes(x = Pruning, y = SNPs, fill = Pruning)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(~Population,nrow = 4,scales = 'free')+
  labs(
       x = "Population",
       y = "SNPs") +
  theme_minimal()

ggsave('BarPlot_pruning.pdf')

to_plot <- table_melt[table_melt$Population %in% c('Founders','A6'),]
ggplot(to_plot, aes(x = Pruning, y = SNPs, fill = Pruning)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(~Population,nrow = 1)+
  labs(
    x = "Population",
    y = "SNPs") +
  theme_minimal()

ggsave('Reduce_BarPlot_pruning.pdf',height = 5)
