res <- read_delim("datasets/WGCNA_data_ClinicalTraits.csv")

library(GGally)

res$sex <- as.character(res$sex)

p <- ggpairs(res,                 # Data frame
             columns = 3:22,        # Columns
             aes(fill = sex))

ggsave(filename = "EDA_prectice_mice.png", p,
       width = 15, height = 15, units = "in", dpi = 300)
