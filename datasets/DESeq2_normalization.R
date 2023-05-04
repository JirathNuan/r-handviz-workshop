library(DESeq2)
library(tidyverse)
library(pcaExplorer)

# 1. Import and concatenate file ====
## 1.1 Load quantification table ====
count.df <- read_delim("datasets/E-GEOD-122495-raw-counts.tsv")

# convert to matrix
count.mat <- as.matrix(count.df[ ,3:8])
rownames(count.mat) <- count.df$`Gene ID`


# 2. Create metadata  ====
meta_df <- read_delim("datasets/E-GEOD-122495-experiment-design.tsv")

# factorize data
meta_df$genotype <- factor(meta_df$genotype,
                           levels = c("Col-0", "pp7l-1"))

# 3. Sample QC with DESeq2 ====
## Generate DESeq object
dds <- DESeqDataSetFromMatrix(round(count.mat), 
                              colData = meta_df, 
                              design = ~ genotype)
dds <- DESeq(dds)

# Create gene dictionary
gene_dict <- select(count.df, `Gene ID`, `Gene Name`)


## Export normalized count
normalized_counts <- counts(dds, normalized=TRUE)

# Rename column name from SRR.... to genotype
colnames(normalized_counts)[colnames(normalized_counts) == "SRR8186097"] <- "pp7l-1_3"
colnames(normalized_counts)[colnames(normalized_counts) == "SRR8186093"] <- "Col-0_2"
colnames(normalized_counts)[colnames(normalized_counts) == "SRR8186095"] <- "pp7l-1_1"
colnames(normalized_counts)[colnames(normalized_counts) == "SRR8186096"] <- "pp7l-1_2"
colnames(normalized_counts)[colnames(normalized_counts) == "SRR8186092"] <- "Col-0_1"
colnames(normalized_counts)[colnames(normalized_counts) == "SRR8186094"] <- "Col-0_3"

# Create data frame
norm_df <- data.frame(normalized_counts) %>% 
  rownames_to_column(var = "Gene ID") %>% 
  left_join(gene_dict, by = "Gene ID")

norm_df <- select(norm_df, `Gene ID`,  `Gene Name`, `Col.0_1`, `Col.0_2`, `Col.0_3`,
         `pp7l.1_1`, `pp7l.1_2`, `pp7l.1_3`)

write_csv(norm_df, "datasets/E-GEOD-122495-normalized_count.csv")

# Get normalized count with collapsed replicate
normCount_collapsed <- collapseReplicates(dds, 
                                          groupby = dds$genotype, 
                                          run = dds$Run)
normCount_collapsed_mat <- counts(normCount_collapsed, normalized=TRUE)


norm_collapsed_df <- data.frame(normCount_collapsed_mat) %>% 
  rownames_to_column(var = "Gene ID") %>% 
  left_join(gene_dict) %>% 
  relocate(`Gene Name`, .after = `Gene ID`)
write_csv(norm_collapsed_df, "datasets/E-GEOD-122495-normalized_count_collapsedRep.csv")




