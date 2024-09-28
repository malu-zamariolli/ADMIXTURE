library(data.table)

#Input files
initial_bim <- fread(snakemake@input[["initial_bim"]], header = FALSE)
final_bim <- fread(snakemake@input[["final_bim"]], header = FALSE)
initial_fam <- fread(snakemake@input[["initial_fam"]], header = FALSE)
final_fam <- fread(snakemake@input[["final_fam"]], header = FALSE)

# print information
text <- paste0("Initial number of individuals: ", nrow(initial_fam), 
               " Initial number of variants: ", nrow(initial_bim),
               " Final number of individuals: ", nrow(final_fam), 
               " Final number of SNPs: ", nrow(final_bim))

# Save
writeLines(text, snakemake@output[["report"]])