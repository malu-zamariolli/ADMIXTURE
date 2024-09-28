
########### Plot admixture per K #################################
library(data.table)
library(dplyr)
library(ggplot2)

###########Input files ###########################

# Obtaining K 
filename <- snakemake@input[["admixture"]]
k <- as.numeric(gsub(".*\\.(\\d+)\\..*", "\\1", filename))
print(paste0("starting K = ", k))

population <- fread(snakemake@input[["population"]], select = c(1, 3))
admixture <- fread(snakemake@input[["admixture"]], col.names = paste0("Cluster", 1:k))
fam <- fread(snakemake@input[["fam"]], select = 2, col.names = "ID")

######### output files ###################
table <- snakemake@output[["table"]]
poptable <- snakemake@output[["poptable"]]
plot <- snakemake@output[["plot"]]


####### Prepare dataframe #############
# Merge df
df_admixture <- cbind(fam, admixture)
df <- merge(df_admixture, population, by = "ID") 

fwrite(df, table,  
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Proportions per groups
result <- df %>%
  group_by(GROUP) %>%
  summarize(across(starts_with("Cluster"), 
                   list(Mean = mean, SD = sd), .names = "{col}_{fn}"))


fwrite(result, poptable,
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

####################### PLotting ###########################################################
# Adapted from https://github.com/mkanai/cancritls/blob/master/admixture/admixture.plot.r

admixture.ordered <- df[order(df$GROUP, decreasing = TRUE), ]
admixture.index <- cbind(index = seq(nrow(admixture.ordered)), admixture.ordered)
admixture.melted <- melt(admixture.index, id = c("index", "ID", "GROUP"), 
                         variable.name = "population")


# Coordinates for annotations
mean_indices <- aggregate(index ~ GROUP, admixture.melted, FUN = mean)
label_pop = mean_indices$GROUP
coord_label <- mean_indices$index
min_indices <- aggregate(index ~ GROUP, admixture.melted, FUN = min)
coord_line <- min_indices$index

p <- ggplot(admixture.melted, aes(x = index, y = value, fill = population)) + 
  geom_bar(stat = "identity", width = 1) +
  labs(fill = paste0("K = ", k)) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  theme_minimal() +
  annotate("text", x = coord_label, y = -0.05, label = label_pop, size = 4) +
  geom_vline(xintercept = coord_line, 
             linetype = "dashed", color = "black")

ggsave(plot, p, height = 12, width = 20)
