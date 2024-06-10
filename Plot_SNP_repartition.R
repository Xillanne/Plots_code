### Script to plot SNPs along chromosomes find by GATK pipelines

setwd("~/Documents/Worms/Plot_GATK/Pop")

All <- read.csv("A6 _snp_positions.csv")
exclude <- read.csv('SNP_2exclude_All.tsv',header = FALSE,sep = '\t',col.names = colnames(All))

library(ggplot2)
library(gridExtra)

All$key <- paste(All$CHROM,All$POS,sep = '_')
exclude$key <- paste(exclude$CHROM,exclude$POS,sep = '_')

All$Filter <- as.factor(ifelse(All$key %in% exclude$key, 'FILTER','PASS'))

All$Plot <- rep(1,nrow(All))

All_filter <- All[All$Filter=='PASS',]
ggplot(All_filter,aes(x=POS, y=Plot,color=Filter))+
  geom_point(size = 0.5)+
  facet_wrap(~ CHROM, ncol = 2, scales = "free_x")+
  theme_bw()+
  ggtitle('SNP repartition for All populations along chromosomes')

nbr <- c()
lenght_chrom <- c()
for(chrom in unique(All_filter$CHROM)){
  nbr <- c(nbr,nrow(All_filter[All_filter$CHROM==chrom,]))
  lenght_chrom <- c(lenght_chrom,All_filter[All_filter$CHROM==chrom,][nrow(All_filter[All_filter$CHROM==chrom,]),]$POS)
}

Summary <- data.frame(Chromosome = unique(All$CHROM),
                      Number_of_SNPs = nbr,
                      Lenght_chromosome = lenght_chrom,
                      Proportion = round(nbr/lenght_chrom*100,5))

g<-ggplot(Summary)+
  geom_point(aes(x=Chromosome,y=Lenght_chromosome,color=Chromosome)) +
  theme_minimal()+
  labs(y='Lenght of chromosome')+
  theme(legend.position = "none", axis.text.x =element_blank(),axis.title.x = element_blank())

g1<-ggplot(Summary)+
  geom_bar(stat = "identity",aes(y = Number_of_SNPs,x=Chromosome,color=Chromosome,fill=Chromosome))+
  geom_text(aes(y = Number_of_SNPs,x=Chromosome,label=paste(Proportion,'%')),vjust=-0.5,size=3)+
  theme_minimal()+
  labs(y='Number of SNPs')+
  theme(legend.position = "none")

plots <- list(g,g1)
pdf('Distribution_SNP_pop.pdf')
grid.arrange(grobs = plots, ncol = 1,heights=c(1,2))
dev.off()
g1
