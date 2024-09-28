library(data.table)

#Input files
final_bim <- fread(snakemake@input[["final_bim"]], header = FALSE)
final_fam <- fread(snakemake@input[["final_fam"]], header = FALSE)

# print information
text <- paste0(" N of individuals after merge: ", nrow(final_fam),
               " N of SNPs after merge: ", nrow(final_bim))

# Save
writeLines(text, snakemake@output[["report"]])