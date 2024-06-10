## Script to create a table which contains all informations about population origin
#To change
setwd("~/Documents/Worms/DataAnalysis/Analyse_summary_tables")

## Reads table where the name of lines are
table_RILs_Tom<-read.delim("summary_reads.RILs.Tom.tsv")
table_CeMee<-read.delim("summary_reads.without_founders.Mapping_RILs.tsv")
table_CeMee_missing<-read.delim("summary_reads.SRR_missing.tsv")
table_00 <- read.delim("summary_reads.tsv")
table_01 <- read.delim("summary_reads_01.tsv")
table_02 <- read.delim("summary_reads_02.without_founders.tsv")

## Set the batch name
liste_table <- list(table_RILs_Tom,table_CeMee,table_CeMee_missing,table_00,table_01,table_02)
batch_list <- c("Recombinant","CeMee","CeMee","NYC_1","NYC_2","NYC_3")


## Initialize Df 
Df_no_concat <- data.frame(Population = character(),
                 Num_generation = character(),
                 Num_line = character(),
                 Batch = character(),
                 stringsAsFactors = FALSE)

Df <- data.frame(Population = character(),
                           Num_generation = character(),
                           Num_line = character(),
                           Batch = character(),
                           stringsAsFactors = FALSE)


for (i in seq_along(liste_table)) {
  # Extract the current table from the list
  table <- liste_table[[i]]
  
  # Extract the corresponding batch name
  batch_name <- batch_list[i]
  Df_no_concat <- extract_line_no_concat(table,Df_no_concat,batch_name)
  Df <- extract_line(table,Df,batch_name)
}

## Sort dataframes
Df <- Df[order(Df$Population), ]
Df_no_concat <- Df_no_concat[order(Df_no_concat$Population), ]

## get the number of replicat
num_rep <- length(rownames(Df[grepl(",",Df$Batch),]))

## Write results in table
write.table(x = Df,file = "Table_composition.tsv",sep = "\t",row.names = FALSE,quote = FALSE)

table(Df_no_concat[c("Population","Batch")]) #Pb on a N2_AN dans CeMee et CB ?
write.table(table(Df_no_concat[c("Population","Batch")]),file = "Table_Batch_X_Pop.tsv",sep = '\t',quote = FALSE)



####### FUNCTIONS #######@


### Function which extract pop, number of gen and line from the name of a line
## if a line is present in two batch it'll change the batch name by putting both of thme
extract_line <- function(table,Df,name_batch){
  Lines <- table[["Line"]]
  # Loop through the elements of Lines
  for (i in Lines) {
    if( !substr(i,1,2) %in% c("A6","CA","GA","GM","GT","LR","SM")){
      pop <- gsub("[0-9].*", "", i) #Get the first letter
      gen <- gsub("^[A-Za-z]+([0-9]+).*", "\\1",i) #Get number after the first letter
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch    
    }
    if (substr(i,1,2)=="A6"){
      pop <- substr(i, 1, 2)
      gen <- substr(i, 3, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch
    }
    if (substr(i,1,2) %in% c("CA","GA","GM","LR","GT")){
      pop <- substr(i, 1, 3)
      gen <- substr(i, 4, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch
    }
    if (substr(i,1,3)=="SMR"){
      pop <- substr(i, 1, 4)
      gen <- NaN
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch 
    }
    if (substr(i,1,2)=="LR"){
      gen <- NaN
    }
    
    #print(paste("pop:", pop, "gen:", gen, "n_line:", n_line, "batch:", batch))
    # Create a new row with the extracted values
    matching_row <- Df$Population == pop & Df$Num_line == n_line & Df$Num_generation == gen
    #print(matching_row)
    if(any(matching_row)){
      Df$Batch[matching_row] <- paste0(Df$Batch[matching_row],paste0(", ",name_batch))
    } else {
      new_row <- data.frame(Population = pop,
                            Num_generation = gen,
                            Num_line = n_line,
                            Batch = batch)
      #print(new_row)
      # Add the new row to the existing data frame
      Df <- rbind(Df, new_row)
    }
  }
  return(Df)
}

### Same function than the previous one but this time even if the line already presents it will put a new row
extract_line_no_concat <- function(table,Df,name_batch){
  Lines <- table[["Line"]]
  # Loop through the elements of Lines
  for (i in Lines) {
    if( !substr(i,1,2) %in% c("A6","CA","GA","GM","GT","LR","SM")){
      pop <- gsub("[0-9].*", "", i) #Get the first letter
      gen <- gsub("^[A-Za-z]+([0-9]+).*", "\\1",i) #Get number after the first letter
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch    
    }
    if (substr(i,1,2)=="A6"){
      pop <- substr(i, 1, 2)
      gen <- substr(i, 3, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch
    }
    if (substr(i,1,2) %in% c("CA","GA","GM","LR","GT")){
      pop <- substr(i, 1, 3)
      gen <- substr(i, 4, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch
    }
    if (substr(i,1,3)=="SMR"){
      pop <- substr(i, 1, 4)
      gen <- NaN
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch 
    }
    if (substr(i,1,2)=="LR"){
      gen <- NaN
    }
    #print(paste("pop:", pop, "gen:", gen, "n_line:", n_line, "batch:", batch))
    new_row <- data.frame(Population = pop,
                          Num_generation = gen,
                          Num_line = n_line,
                          Batch = batch)
    #print(new_row)
    # Add the new row to the existing data frame
    Df <- rbind(Df, new_row)
  }
  return(Df)
}

