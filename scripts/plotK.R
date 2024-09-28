########### Plotting K and CV errors from admixture ###############
library(data.table)
library(ggplot2)

#################################################################

# input data
cv <- fread(snakemake@input[["cv"]], col.names = c("K", "CVerror"))

cv$K <- gsub("\\(K=", " ", cv$K)
cv$K <- gsub("):", "", cv$K)
cv$K <- as.numeric(cv$K)

cv <- cv[order(cv$K), ]

kmin <- as.numeric(cv[1, "K"])
kmax <- as.numeric(cv[nrow(cv), "K"])

# Plot
p <- ggplot(cv, aes(x = K, y = CVerror)) +
  geom_line() +
  scale_x_continuous(breaks = c(kmin:kmax), labels = c(kmin:kmax)) +
  ylab("Cross-validation error") +
  theme_minimal()

ggsave(snakemake@output[["plot"]], p, height = 8, width = 8)