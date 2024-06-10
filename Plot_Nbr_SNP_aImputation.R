### Script pour comparer le nombre de SNPs par population après imputation

setwd("~/Documents/Worms/Plot_GATK/Pop")



library(tidyr)
library(ggplot2)

# Charger les fichiers CSV
files <- list.files(pattern = ".*.imputed _snp_positions.csv")
files_name <- gsub('founders','F',sub("(.*).imputed _snp_positions.csv",'\\1',files))
dfs <- lapply(files, read.csv)

nbr <- c()
for(i in 1:length(files)){
  nbr <- c(nbr,nrow(dfs[[i]]))
}

To_plot <- data.frame(nbr = nbr,
                      pop = files_name)

ggplot(To_plot)+
  geom_bar(stat = "identity",aes(y = nbr,x=pop,color=pop,fill=pop))+
  geom_text(aes(y = nbr,x=pop,label=nbr),vjust=-0.5)+
  theme_minimal()+
  labs(y='Number of SNPs')+
  theme(legend.position = "none")

ggsave('Number_SNPs_aImputation.pdf')

library(combinat)
library(dplyr)

comb <- combn(1:length(files),2)
comb_name <- combn(files_name,2)

nbr_common <- c()
name <- c()

for(i in 1:ncol(comb)){
  df1 <- dfs[[comb[,i][1]]]
  df2 <- dfs[[comb[,i][2]]]
  name <-c(name,paste(comb_name[,i][1],comb_name[,i][2],sep = ','))
  if(nrow(df1)==0 || nrow(df2)==0){
    nbr_common<-c(nbr_common,0)
  }else{
  nbr_common <- c(nbr_common,nrow(intersect(df1,df2)))}
}


data <- data.frame(Comparison = name,
                   Nbr_SNPs = nbr_common,
                   nbr_string = sprintf("%.2e",round(nbr_common)))

# On retir les valeurs à 0
data <- data[data$Nbr_SNPs!=0,]

ggplot(data, aes(x = Comparison, y = Nbr_SNPs, fill = Comparison)) +
  geom_bar(stat = "identity") +
  geom_point()+
  labs(x = "Traits", y = "Number of SNPs", title = "Number of SNPs Shared between Pop (EEV removed)") +
  theme_minimal()+
  geom_text(aes(x = Comparison, y = Nbr_SNPs,label=Nbr_SNPs),vjust=-0.5,size=1.75)+
  scale_x_discrete(labels = function(x) gsub(",", "\n", x))+
  theme(legend.position = "none", axis.text.x = element_text(size = 7))

ggsave('Shared_SNPs_pop.pdf')


# Séparer les noms de populations pour créer une matrice

data_split <- separate(data, Comparison, into = c("Pop1", "Pop2"), sep = ",")
data_split$Pop1 <- factor(data_split$Pop1, levels = c('F','A6','CA','GA','GM','GT','LR','SMR'))
data_split$Pop2 <- factor(data_split$Pop2, levels = c('F','A6','CA','GA','GM','GT','LR','SMR'))
data_split$nbr_string <- NULL

# Ajouter les valeurs diagonales
diagonal_data <- data.frame(Pop1 = files_name, Pop2 = files_name, Nbr_SNPs = nbr)
diagonal_data <- diagonal_data[diagonal_data$Pop1 != 'EEV', ]

heatmap_data <- rbind(data_split, diagonal_data)
heatmap_data[heatmap_data$Pop2=='F' & heatmap_data$Pop1=='A6', ][,c('Pop1','Pop2')] <- c('F','A6')
heatmap_data[heatmap_data$Pop2=='F' & heatmap_data$Pop1=='CA', ][,c('Pop1','Pop2')] <- c('F','CA')

# Créer la heatmap
ggplot(heatmap_data, aes(x = Pop1, y = Pop2, fill = Nbr_SNPs)) +
  geom_tile() +
  scale_fill_gradient(low = 'yellow', high = "red") +
  geom_text(aes(label = Nbr_SNPs), color = "black", size = 2) +
  theme_minimal() +
  labs(x = "Population 1", y = "Population 2", fill = "Number of SNPs", title = "Heatmap of Shared SNPs Between Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('Heatmap_shared_SNPs.pdf')
