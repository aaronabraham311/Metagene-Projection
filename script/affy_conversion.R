library(bioMart)
library(tidyr)
library(mouse4302.db) # Annotation data for affymetrix mouse array 430 2.0

#Reads in data
df <- read.delim("data/control2.txt", sep = "\t")
ortho_human_genes <- read.delim("data/mart_export.txt")

df_expression <- df
colnames(df_expression) <- c("probe_id", "expression")
colnames(ortho_human_genes) <- c("ensembl_id", "transcript_id", "human_id")

# Maps Affymetrix data into Ensembl ID
df = toTable(mouse4302ENSEMBL)
new_df <- merge(df, df_expression, by = "probe_id") # Combining expression data based on common probe_id
write.table(new_df, "data/controlNew_2.txt", sep = "\t", row.names = F)