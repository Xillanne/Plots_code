### Script for analysing PSC results from output of bcftools stats

# First PSC need to be grep from output of bcftools stats with the following cmd
#grep "^PSC" output.vchk > output_PSC.tsv

# Name working directory 
setwd("~/Documents/Worms/VariantCalling/Stats")
# Name of the grepped file
file <- "output_AllLines.bF_PSC.vchk" 
# Name of output
name <- "AllLines.bF.sampleRemoved"
# Cutoff : extract lines with values >
T_singletons <- 150
T_het <- 0.085

library(ggplot2)
library("ggsci")

# Convert file in dataframe format
df <- read.delim(file,header = FALSE,)
to_removed <- read.csv('/Users/alix/Documents/Worms/VariantCalling/Isotypes/Removed_line.csv')

df <- df[!(df$V3 %in% to_removed$Line),]

# Set columns name
colnames(df) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions",	"nIndels","average_depth","nSingleton","nHapRef","nHapAlt","nMissing")
# Keep only columns we are intersseted in
Df <- df[c("sample","nRefHom","nNonRefHom","nHets","nSingleton","nMissing")]

# Remove old dataframe
rm(df)

# Distinc different pop
populations<-c("A6","CA1","CA2","CA3","EEV","GA1","GA2","GA4","GM1","GM3","GT1","GT2","LR1","LR3","SMR2","SMR4")
Df$population <- ifelse(grepl("A6", Df$sample), "A6",
                 ifelse(grepl("CA1", Df$sample), "CA1",
                        ifelse(grepl("CA2", Df$sample), "CA2",
                               ifelse(grepl("CA3", Df$sample), "CA3",
                                      ifelse(grepl("EEV", Df$sample), "EEV",
                                             ifelse(grepl("GA1", Df$sample), "GA1",
                                                    ifelse(grepl("GA2", Df$sample), "GA2",
                                                           ifelse(grepl("GA4", Df$sample), "GA4",
                                                                  ifelse(grepl("GM1", Df$sample), "GM1",
                                                                         ifelse(grepl("GM3", Df$sample), "GM3",
                                                                                ifelse(grepl("GT1", Df$sample), "GT1",
                                                                                       ifelse(grepl("GT2", Df$sample), "GT2",
                                                                                              ifelse(grepl("LR1", Df$sample), "LR1",
                                                                                                     ifelse(grepl("LR3", Df$sample), "LR3",
                                                                                                            ifelse(grepl("SMR2", Df$sample), "SMR2",
                                                                                                                   ifelse(grepl("SMR4", Df$sample), "SMR4", NA))))))))))))))))

if(all(is.na(Df$population))){Df$population<-Df$sample}
# Plot Number of SNP per sample
Df$nSNP <- Df$nNonRefHom + Df$nHets #nRefHom is not a SNP because not polymorphism
Df$sampleN <- seq(0,length(Df$sample)-1)
g<-ggplot()+
  geom_point(data=Df,aes(x=sampleN,y=nSNP, color=population))+
  ylab('Number of SNP')+
  ggtitle(gsub("XXX",name,"Number of SNP per sample XXX"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 


# Plot heterozygotie rate per line
Df$rHetero <- Df$nHets/(Df$nRefHom+Df$nNonRefHom)
h<-ggplot()+
  geom_point(data=Df,aes(x=sampleN,y=rHetero, color=population))+
  ylab('Heterozygotie fraction (nHet/nHom)')+
  ggtitle(gsub("XXX",name,"Heterozygotie rate per sample XXX"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 



# Plot Number of singleton per line
s<-ggplot()+
  geom_point(data=Df,aes(x=sampleN,y=nSingleton, color=population))+
  ylab('Number of singletons')+
  theme_bw()+
  ggtitle(gsub("XXX",name,"Number of singletons per sample XXX"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

# Plot Number of missing value per line
m<-ggplot()+
  geom_point(data=Df,aes(x=sampleN,y=nMissing, color=population))+
  ylab('Number of missing')+
  theme_bw()+
  ggtitle(gsub("XXX",name,"Number of missing values per sample XXX"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

# Save plot as PDF
pdf_filename <- gsub("XXX",name,"Plots_XXX.pdf")

# Open PDF device
pdf(pdf_filename) #, width = 8, height = 4

# Print each plot to the PDF
print(g)
print(h)
print(s)
print(m)

# Close PDF device
dev.off()

# Get Lines with singletons and hetero > at some cutoffs
Lines_singletons <- Df[Df$nSingleton>T_singletons,]$sample
Lines_het <- Df[Df$rHetero>T_het,]$sample
# remove Na values
Lines_singletons <- na.omit(Lines_singletons)
Lines_het <- na.omit(Lines_het)
#order
Lines_singletons <- Lines_singletons[order(Lines_singletons)]
Lines_het <- Lines_het[order(Lines_het)]

write.csv(data.frame(Line=Lines_het,Issue=rep('High_heterozygotity',length(Lines_het))),'Removed_line_hetero.csv',row.names = FALSE,quote = FALSE)

#Get lines with none SNPs 
Lines_none <- Df[Df$nSNP==0,]
Lines_none_0 <- Lines_none[Lines_none$nRefHom==0,][c('sample','nMissing')]
Lines_none_idem <- Lines_none[Lines_none$nRefHom!=0,][c('sample','nMissing','nRefHom')]


