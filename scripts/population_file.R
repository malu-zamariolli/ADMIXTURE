# Create file with populations from all samples in the final file

library(data.table)

# Files
fam <- fread(snakemake@input[["fam"]], header = FALSE, select = c(2), col.names = "ID")
fambr <- fread(snakemake@input[["fambr"]], header = FALSE, select = c(2), col.names = "ID")
pop1k <- fread(snakemake@input[["pop1k"]], header = TRUE, select = c(1:3))

output <- snakemake@output[["file"]]

# Create columns for episono cohort
fambr$POP <- c("BRSP")
fambr$GROUP <- c("ADMIXED")

# Rbind episono and 1K population file
df <- rbind(fambr, pop1k)

# Extract samples in fam from df
df_final <- df[df$ID %in% fam$ID, ]

print(paste0("Dimension of merged fam file: ", dim(fam)))
print(paste0("Dimension of population file: ", dim(df_final)))

fwrite(df_final, output, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


