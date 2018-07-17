# Reading libraries
library(tidyr)
library(caret)
library(data.table)

# Reading data
df <- read.delim("TPM_table.txt", sep = "\t", row.names = 1)

# Removing rows with 0 values (low-expression genes)
df$X <- NULL
row_sub = apply (df[, 1:8], 1, function(row) all(row != 0))
newDF <- df[row_sub,]

# Removing low variance features (i.e. all values are between 0.5 standard deviations)
row_sub = apply(newDF, 1, function(row) all(var(row) > 0.5))
newDF <- newDF[row_sub,]

# Putting row names in column and writing dataset
setDT(newDF, keep.rownames = TRUE)
colnames(newDF)[colnames(newDF) == "rn"] <- "Genes"
write.table(newDF, "model_data.txt", sep = "\t")
