library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)

# loading in the DESeq object
dds <- readRDS("dds.rds")

### PCA ###

# Regularized log transformations for PCA
rld <- rlog(dds)

# Generating PCA data
pcaData <- plotPCA(rld, intgroup=c("Timepoint"), returnData = TRUE)
pcaData

# Extracting the variation percents from the PCA object
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA Plot
pcaData %>% 
  ggplot(aes(x = PC1, y = PC2, colour = Timepoint)) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(-12,12)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_test() +
  theme(
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(linewidth = 2)) +
  labs(colour = "Timepoint")

### HEATMAP ###

# Removing ncRNA genes from dds object
ncrna_genes <- grepl("NCRNA", rownames(dds), ignore.case = F)

dds_filtered <- dds[!ncrna_genes, ]

# Removing SPOM_ prefix from gene names
rownames(dds_filtered) <- gsub("SPOM_", "", rownames(dds_filtered))

# Perform log transformation on filtered count data
rld_1 <- rlog(dds_filtered)
rld_mat <- assay(rld_1)

# Selecting top 20 most variable genes (can be altered to generate top 500)
select <- head(order(rowVars(rld_mat), decreasing = TRUE), 20)
selected_genes <- assay(rld_1[select,])
selected_genes_names <- rownames(selected_genes)

# Loading in gene id names from PomBase
pom_genes <- read_tsv('data/gene_IDs_names.tsv', skip = 1)

# Converting gene IDs to gene names
filtered_pom_genes <- pom_genes[tolower(pom_genes$gene_systematic_id) %in% tolower(selected_genes_names), ]

filtered_pom_genes_vec <- ifelse(is.na(filtered_pom_genes$gene_name), 
                                 filtered_pom_genes$gene_systematic_id, 
                                 filtered_pom_genes$gene_name)

rownames(selected_genes) <- filtered_pom_genes_vec

# Data frame of timepoint names for heatmap
df <- as.data.frame(colData(dds_filtered)['Timepoint'])

# Getting gene names for plot
rownames(df) <- colnames(rld_1)

# Generating heatmap
pheatmap(selected_genes, scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=F, annotation_col=df)



