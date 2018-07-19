# Reading libraries
library(tidyr)
library(plyr)
library(caret)
library(data.table)
library(bioMart)
library(hgu133a.db)

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
write.table(all_df, "data/uncleaned-human-data.txt", sep = '\t')

#Converting Affymetrix data into ENSEMBL gene ids
all_df[,1] <- gsub('.*:',"", all_df[,1])
colnames(all_df)[1] <- "probe_id"
ensemblKey = toTable(hgu133aENSEMBL)

all_df <- merge(all_df, ensemblKey, by = "probe_id")
all_df$probe_id <- all_df$ensembl_id
all_df$ensembl_id <- NULL
colnames(all_df)[1] <- "Genes"

# Averaging out duplicated genes
df_new <- data.table(all_df)
df_new <- df_new[,list(ARS_3 = mean(ARS_3), ARS_4 = mean(ARS_4), ERS_1 = mean(ERS_1),
                       ERS_2 = mean(ERS_2), CBS_1 = mean(CBS_1), CBS_2 = mean(CBS_2), 
                       ES_1 = mean(ES_1), ES_2 = mean(ES_2), ES_3 = mean(ES_3), FMT_1 = mean(FMT_1), 
                       FMT_2 = mean(FMT_2), MSS_1 = mean(MSS_1), MSS_2 = mean(MSS_2), MSS_3 = mean(MSS_3), 
                       MSS_4 = mean(MSS_4), LMS_1 = mean(LMS_1), LMS_2 = mean(LMS_2), LMS_3 = mean(LMS_3),
                       LMS_4 = mean(LMS_4), MLS_1 = mean(MLS_1), MLS_2 = mean(MLS_2), MLS_3 = mean(MLS_3), 
                       MLS_4 = mean(MLS_4), OS_1 = mean(OS_1), OS_2 = mean(OS_2), OS_3 = mean(OS_3), 
                       OS_4 = mean(OS_4), OS_5 = mean(OS_5), OS_6 = mean(OS_6), ARS_1 = mean(ARS_1),
                       ERS_3 = mean(ERS_3), CBS_3 = mean(CBS_3), CBS_4 = mean(CBS_4), ES_4 = mean(ES_4),
                       ES_5 = mean(ES_5), FMT_3 = mean(FMT_3), FMT_4 = mean(FMT_4), MSS_5 = mean(MSS_5), 
                       MSS_6 = mean(MSS_6), MSS_7 = mean(MSS_7), MSS_8 = mean(MSS_8), LMS_5 = mean(LMS_5), 
                       LMS_6 = mean(LMS_6), LMS_7 = mean(LMS_8), MLS_5 = mean(MLS_5), MLS_6 = mean(MLS_6),
                       MLS_7 = mean(MLS_7), OS_7 = mean(OS_7), OS_8 = mean(OS_8), OS_9 = mean(OS_9), 
                       OS_10 = mean(OS_10), OS_11 = mean(OS_11)), list(Genes)]
all_df <- df_new
write.table(all_df, "data/human_gene_data.txt", sep = "\t")

# Combining human genes and mouse genes
h_df <- read.delim("data/human_gene_data.txt", sep = "\t")
m_df <- read.delim("data/mouse_gene_data.txt", sep = "\t")

h_df <- h_df[h_df$Genes %in% m_df$Genes, ]
m_df <- m_df[m_df$Genes %in% h_df$Genes, ] 

combined_df <- merge(h_df, m_df, by = "Genes")

write.table(combined_df, "data/combined_data.txt", sep = "\t")

# Removing rows with 0 values (low-expression genes)
df$X <- NULL
row_sub = apply (df[, 2:12], 1, function(row) all(row != 0))
newDF <- df[row_sub,]

# Removing low variance features (i.e. all values are between 0.5 standard deviations)
row_sub = apply(newDF, 1, function(row) all(var(row[2:12]) > 0))
newDF <- newDF[row_sub,]

# Putting row names in column and writing dataset
write.table(newDF, "model_data.txt", sep = "\t")
