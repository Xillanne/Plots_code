## SCRIPT FOR ANALYSING TABLE OF READ SUMMARY FOR SEQUENCING ##

# Directory where the tables are
setwd("~/Documents/Worms/DataAnalysis/Analyse_summary_tables")

# Librairies
library(ggplot2)


read_table <- function(filename) {
  # Read the table
  table <- read.table(filename, header = TRUE)
  
  # Calculate percentages
  table$percentage_after_fastp <- table$Reads_after_fastp * 100 / table$Reads
  table$percentage_final <- table$Reads_final * 100 / table$Reads
  
  return(table)
}

table_RILs_Tom<-read_table("summary_reads.RILs.Tom.tsv")
table_CeMee<-read_table("summary_reads.without_founders.Mapping_RILs.tsv")
table_CeMee_founders<-read_table("summary_reads.founders.Mapping_RILs.tsv")
table_00 <- read_table("summary_reads.tsv")
table_01 <- read_table("summary_reads_01.tsv")
table_02 <- read_table("summary_reads_02.without_founders.tsv")
table_founders_NYC <- read_table("summary_reads_02.founders.tsv")
table_CeMee_missing <- read_table("summary_reads.SRR_missing.tsv")

NYC <- data.frame(rbind(table_00,
                        table_01,
                        table_02))

CeMee <- data.frame(rbind(table_CeMee,
                        table_CeMee_missing))

CeMee$Line <- paste0(CeMee$Line,'CeMee')

All <- data.frame(rbind(NYC,
                        CeMee,
                        table_RILs_Tom))

write.table(All,'AllLines.summary_table.tsv',sep = '\t',quote = FALSE,row.names = FALSE)


# Plot the points with specified color
## 
ggplot()+
  geom_point(data=table_RILs_Tom, aes(x = percentage_final, y = percentage_after_fastp, color = "Rils_Tom")) +
  geom_point(data=table_CeMee, aes(x = percentage_final, y = percentage_after_fastp, color = "CeMee"))+
  geom_point(data=table_CeMee_missing, aes(x = percentage_final, y = percentage_after_fastp, color = "CeMee"))+
  geom_point(data=table_CeMee_founders, aes(x = percentage_final, y = percentage_after_fastp, color = "Founders_CeMee"))+#ylim(75,100)+
  geom_point(data=table_00, aes(x = percentage_final, y = percentage_after_fastp, color = "NYC_batch1"))+
  geom_point(data=table_01, aes(x = percentage_final, y = percentage_after_fastp, color = "NYC_batch2"))+
  geom_point(data=table_02, aes(x = percentage_final, y = percentage_after_fastp, color = "NYC_batch3"))+
  geom_point(data=table_founders_NYC, aes(x = percentage_final, y = percentage_after_fastp, color = "Founders_NYC"))+
  #geom_point(data=table_founders_NYC, aes(x = percentage_final, y = percentage_after_fastp, color = "Founders_NYC"))+
  scale_color_manual(values = c(Rils_Tom = "#17d843", CeMee = "#fc0000", Founders_CeMee = "#ffc171",NYC_batch1="#252ff8",NYC_batch2="#2585f8",NYC_batch3="#25b5f8",Founders_NYC="#991bf0")) +
  labs(color = "Group") +  # Customize legend title
  ggtitle('Representation of line by their percentage of reads keep for each sequencing groups')+
  geom_vline(xintercept = 50, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 78, linetype = "dashed", color = "red")+
  theme_bw() 

ggplot()+
  geom_point(data=table_RILs_Tom, aes(x = percentage_final, y = Coverage, color = "Rils_Tom")) +
  geom_point(data=table_CeMee, aes(x = percentage_final, y = Coverage, color = "CeMee"))+
  geom_point(data=table_CeMee_missing, aes(x = percentage_final, y = Coverage, color = "CeMee"))+
  geom_point(data=table_CeMee_founders, aes(x = percentage_final, y = Coverage, color = "Founders_CeMee"))+#ylim(75,100)+
  geom_point(data=table_00, aes(x = percentage_final, y = Coverage, color = "NYC_batch1"))+
  geom_point(data=table_01, aes(x = percentage_final, y = Coverage, color = "NYC_batch2"))+
  geom_point(data=table_02, aes(x = percentage_final, y = Coverage, color = "NYC_batch3"))+
  geom_point(data=table_founders_NYC, aes(x = percentage_final, y = Coverage, color = "Founders_NYC"))+
  #geom_point(data=table_founders_NYC, aes(x = percentage_final, y = Coverage, color = "Founders_NYC"))+
  scale_color_manual(values = c(Rils_Tom = "#17d843", CeMee = "#fc0000", Founders_CeMee = "#ffc171",NYC_batch1="#252ff8",NYC_batch2="#2585f8",NYC_batch3="#25b5f8",Founders_NYC="#991bf0")) +
  labs(color = "Group") +  # Customize legend title
  ggtitle('Representation of lines by their coverage in fonction of the percentage of reads kept for each sequencing groups')+
  geom_hline(yintercept = 5, linetype = "dashed", color = "red")+
  #ylim(0,50)+
  theme_bw() 

ggplot()+
  geom_point(data=table_RILs_Tom, aes(x = percentage_final, y = Reads, color = "Rils_Tom")) +
  geom_point(data=table_CeMee, aes(x = percentage_final, y = Reads, color = "CeMee"))+
  geom_point(data=table_CeMee_missing, aes(x = percentage_final, y = Reads, color = "CeMee"))+
  geom_point(data=table_CeMee_founders, aes(x = percentage_final, y = Reads, color = "Founders_CeMee"))+#ylim(75,100)+
  geom_point(data=table_00, aes(x = percentage_final, y = Reads, color = "NYC_batch1"))+
  geom_point(data=table_01, aes(x = percentage_final, y = Reads, color = "NYC_batch2"))+
  geom_point(data=table_02, aes(x = percentage_final, y = Reads, color = "NYC_batch3"))+
  geom_point(data=table_founders_NYC, aes(x = percentage_final, y = Reads, color = "Founders_NYC"))+
  #geom_point(data=table_founders_NYC, aes(x = percentage_final, y = Reads, color = "Founders_NYC"))+
  scale_color_manual(values = c(Rils_Tom = "#17d843", CeMee = "#fc0000", Founders_CeMee = "#ffc171",NYC_batch1="#252ff8",NYC_batch2="#2585f8",NYC_batch3="#25b5f8",Founders_NYC="#991bf0")) +
  labs(color = "Group") +  # Customize legend title
  geom_vline(xintercept = 50, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 78, linetype = "dashed", color = "red")+
  theme_bw() 


 # Zoom with the specified threshorld
# Set Threshold and select point to plot
x=50
y=100

zoom_table <- function(table,x,y) {
  table_zoom <- table[table$percentage_after_fastp<y&table$percentage_final<x,]
  return(table_zoom)
}
  

table_00_zoom <- zoom_table(table_00,x,y)
table_01_zoom <- zoom_table(table_01,x,y)
table_02_zoom <- zoom_table(table_02,x,y)
table_CeMee_zoom <- zoom_table(table_CeMee,x,y)
table_CeMee_zoom_missing <- zoom_table(table_CeMee_missing,x,y)
table_CeMee_founders_zoom <- zoom_table(table_CeMee_founders,x,y)
table_RILs_Tom_zoom <- zoom_table(table_RILs_Tom,x,y)
table_founders_NYC_zoom <- zoom_table(table_founders_NYC,x,y)

ggplot()+
  geom_point(data=table_RILs_Tom_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "Rils_Tom")) +
  geom_point(data=table_CeMee_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "CeMee"))+
  geom_point(data=table_CeMee_zoom_missing, aes(x = percentage_final, y = percentage_after_fastp, color = "CeMee"))+
  geom_point(data=table_CeMee_founders_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "Founders_CeMee"))+#ylim(75,100)+
  geom_point(data=table_00_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "NYC_batch1"))+
  geom_point(data=table_01_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "NYC_batch2"))+
  geom_point(data=table_02_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "NYC_batch3"))+
  geom_point(data=table_founders_NYC_zoom, aes(x = percentage_final, y = percentage_after_fastp, color = "Founders_NYC"))+
  #geom_point(data=table_founders_NYC, aes(x = percentage_final, y = percentage_after_fastp, color = "Founders_NYC"))+
  scale_color_manual(values = c(Rils_Tom = "#17d843", CeMee = "#fc0000", Founders_CeMee = "#ffc171",NYC_batch1="#252ff8",NYC_batch2="#2585f8",NYC_batch3="#25b5f8",Founders_NYC="#991bf0")) +
  labs(color = "Group") +  # Customize legend title
  #ggtitle('Zoom of Representation of line by their percentage of reads keep for each sequencing groups')+
  #ylim(0,2.5)+xlim(0,0.2)+
  xlab("percentage_after_mapping")+
  theme_bw() 


ggplot()+
  geom_point(data=table_RILs_Tom_zoom, aes(x = percentage_final, y = Coverage, color = "Rils_Tom")) +
  geom_point(data=table_CeMee_zoom, aes(x = percentage_final, y = Coverage, color = "CeMee"))+
  geom_point(data=table_CeMee_zoom_missing, aes(x = percentage_final, y = Coverage, color = "CeMee"))+
  geom_point(data=table_CeMee_founders_zoom, aes(x = percentage_final, y = Coverage, color = "Founders_CeMee"))+#ylim(75,100)+
  geom_point(data=table_00_zoom, aes(x = percentage_final, y = Coverage, color = "NYC_batch1"))+
  geom_point(data=table_01_zoom, aes(x = percentage_final, y = Coverage, color = "NYC_batch2"))+
  geom_point(data=table_02_zoom, aes(x = percentage_final, y = Coverage, color = "NYC_batch3"))+
  geom_point(data=table_founders_NYC_zoom, aes(x = percentage_final, y = Coverage, color = "Founders_NYC"))+
  #geom_point(data=table_founders_NYC, aes(x = percentage_final, y = Coverage, color = "Founders_NYC"))+
  scale_color_manual(values = c(Rils_Tom = "#17d843", CeMee = "#fc0000", Founders_CeMee = "#ffc171",NYC_batch1="#252ff8",NYC_batch2="#2585f8",NYC_batch3="#25b5f8",Founders_NYC="#991bf0")) +
  labs(color = "Group") +  # Customize legend title
  ggtitle('Zoom of Representation of line by their percentage of reads keep for each sequencing groups')+
  #ylim(y,100)+xlim(0,x)+
  theme_bw() 





line_2delete <- function(table) {
  table_d <- table[table$percentage_after_fastp > 78 & table$percentage_final ==0 & table$Coverage ==0, ]
  return(table_d$Line)
}

list_00_delete <- line_2delete(table_00)
list_01_delete <- line_2delete(table_01)
list_02_delete <- line_2delete(table_02)
list_CeMee_delete <- line_2delete(table_CeMee)
list_CeMee_missing_delete <- line_2delete(table_CeMee_missing)
list_CeMee_founders_delete <- line_2delete(table_CeMee_founders)
list_RILs_Tom_delete <- line_2delete(table_RILs_Tom)
list_founders_NYC_delete <- line_2delete(table_founders_NYC)
  