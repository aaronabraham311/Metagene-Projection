# Reading libraries
library(tidyr)
library(caret)
library(data.table)

# Reading data and cleaning
df <- read.delim("data/TPM_table.txt", sep = "\t", row.names = 1)
setDT(df, keep.rownames = TRUE)
colnames(df)[colnames(df) == "rn"] <- "Genes"

# Change transcripts into genes for original data
mart <- read.delim("data/mart_export.txt", sep = '\t')
mart <- lapply(mart, as.character)
df$Genes <- mapvalues(df$Genes, mart$Transcript.stable.ID, mart$Gene.stable.ID)

# Deleting rows that are still transcripts (i.e. they are not orthologous to human genes)
df <- df[-grep("ENSMUST",df$Genes),]

#Averaging out values of duplicate rows
df_new <- data.table(df)
df_new <- df_new[,list(G3M19LEG = mean(G3M19LEG), G3M5CM = mean(G3M5CM), G3M5LEG = mean(G3M5LEG),
                       RMS1 = mean(RMS1), RMS3 = mean(RMS3), RMS_3 = mean(RMS_3), T_LEFT = mean(T_LEFT),
                       T_LEG = mean(T_LEG)), list(Genes)]
df <- df_new

#Reading in controls
control_1 <- read.delim("data/controlNew_1.txt", sep = "\t")
control_2 <- read.delim("data/controlNew_2.txt", sep = "\t")
control_3 <- read.delim("data/controlNew_3.txt", sep = "\t")

#Renaming columns so it is easier to remove later
control_1$probe_id <- NULL
control_2$probe_id <- NULL
control_3$probe_id <- NULL

colnames(control_1) <- c("Genes", "expression_1")
colnames(control_2) <- c("Genes", "expression_2")
colnames(control_3) <- c("Genes", "expression_3")

# Averaging out duplicates
df_new <- data.table(control_1)
df_new <- df_new[,list(expression_1= mean(expression_1)), list(Genes)]
control_1<- df_new

df_new <- data.table(control_2)
df_new <- df_new[,list(expression_2= mean(expression_2)), list(Genes)]
control_2 <- df_new

df_new <- data.table(control_3)
df_new <- df_new[,list(expression_3= mean(expression_3)), list(Genes)]
control_3<- df_new

# Combining controls
df <- df[df$Genes %in% control_1$Genes]
df <- merge(df, control_1, by = "Genes")
df <- merge(df, control_2, by = "Genes")
df <- merge(df, control_3, by = "Genes")

# Removing rows with 0 values (low-expression genes)
df$X <- NULL
row_sub = apply (df[, 2:12], 1, function(row) all(row != 0))
newDF <- df[row_sub,]

# Removing low variance features (i.e. all values are between 0.5 standard deviations)
row_sub = apply(newDF, 1, function(row) all(var(row[2:12]) > 0))
newDF <- newDF[row_sub,]

# Putting row names in column and writing dataset
write.table(newDF, "model_data.txt", sep = "\t")
