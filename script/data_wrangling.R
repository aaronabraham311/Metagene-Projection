# Reading libraries
library(tidyr)
library(plyr)
library(caret)
library(data.table)

# Reading data and cleaning
df <- read.delim("data/TPM_table.txt", sep = "\t", row.names = 1)
setDT(df, keep.rownames = TRUE)
colnames(df)[colnames(df) == "rn"] <- "Genes"

# Change transcripts into genes for original data
mart <- read.delim("data/mart_export.txt", sep = '\t')
mart$Gene.stable.ID <- as.character(mart$Gene.stable.ID)
mart$Transcript.stable.ID <- as.character(mart$Transcript.stable.ID)
mart$Human.gene.stable.ID <- as.character(mart$Human.gene.stable.ID)

df$Genes <- mapvalues(df$Genes, mart$Transcript.stable.ID, mart$Human.gene.stable.ID)

# Deleting rows that are still transcripts (i.e. they are not orthologous to human genes)
df <- df[-grep("ENSMUST",df$Genes),]

#Averaging out values of duplicate rows
df_new <- data.table(df)
df_new <- df_new[,list(G3M19LEG = mean(G3M19LEG), G3M5CM = mean(G3M5CM), G3M5LEG = mean(G3M5LEG),
                       RMS1 = mean(RMS1), RMS3 = mean(RMS3), RMS_3 = mean(RMS_3), T_LEFT = mean(T_LEFT),
                       T_LEG = mean(T_LEG)), list(Genes)]
df <- df_new
write.table(df, "mouse_gene_data.txt", sep = "\t")


#Reading human gene data
model_df <- read.delim("data/model_data.txt", sep = "\t")
test_df <- read.delim("data/test_data.txt", sep = "\t")

i <- sapply(model_df, is.factor)
model_df[i] <- lapply(model_df[i], as.character)
i <- sapply(test_df, is.factor)
test_df[i] <- lapply(test_df[i], as.character)

colnames(model_df) <- model_df[2,]
model_df <- model_df[-c(1,2),]
colnames(test_df) <- test_df[2,]
test_df <- test_df[-c(1,2),]

model_df[-c(1,2)] <- lapply(model_df[-c(1,2)], as.numeric)
test_df[-c(1,2)] <- lapply(test_df[-c(1,2)], as.numeric)
test_df <- test_df[,-c(2)]

all_df <- merge(model_df, test_df, by = "Name")






# Removing rows with 0 values (low-expression genes)
df$X <- NULL
row_sub = apply (df[, 2:12], 1, function(row) all(row != 0))
newDF <- df[row_sub,]

# Removing low variance features (i.e. all values are between 0.5 standard deviations)
row_sub = apply(newDF, 1, function(row) all(var(row[2:12]) > 0))
newDF <- newDF[row_sub,]

# Putting row names in column and writing dataset
write.table(newDF, "model_data.txt", sep = "\t")
