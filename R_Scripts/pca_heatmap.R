library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)

# Load files into R

dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix,
                                     colData = metadata,
                                     design = ~Timepoint)

dds_matrix

dds_matrix$condition <- relevel(dds_matrix$Timepoint, ref = "0_hour")

# Check the levels of dds_matrix$condition
levels(dds_matrix$Timepoint)
dds <- DESeq(dds_matrix)
saveRDS(dds, "dds.rds")

rld <- rlog(dds)

# Generate a PCA plot with DESeq2's plotPCA function
pcaData <- plotPCA(rld, intgroup=c("Timepoint"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

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


ncrna_genes <- grepl("NCRNA", rownames(dds), ignore.case = F)

# Remove these genes from the dataset
dds_filtered <- dds[!ncrna_genes, ]

rownames(dds_filtered) <- gsub("SPOM_", "", rownames(dds_filtered))

# Perform log transformation on our count data
rld_1 <- rlog(dds_filtered)
# vsd <- vst(dds_filtered)
rld_mat <- assay(rld_1)

select <- head(order(rowVars(rld_mat), decreasing = TRUE), 20)
selected_genes <- assay(rld_1[select,])
selected_genes_names <- rownames(selected_genes)

pom_genes <- read_tsv('data/gene_IDs_names.tsv', skip = 1)

filtered_pom_genes <- pom_genes[tolower(pom_genes$gene_systematic_id) %in% tolower(selected_genes_names), ]

filtered_pom_genes_vec <- ifelse(is.na(filtered_pom_genes$gene_name), 
                                 filtered_pom_genes$gene_systematic_id, 
                                 filtered_pom_genes$gene_name)

rownames(selected_genes) <- filtered_pom_genes_vec

df <- as.data.frame(colData(dds_filtered)['Timepoint'])

rownames(df) <- colnames(rld_1)

pheatmap(selected_genes, scale = "row", cluster_rows=T, show_rownames=T,
         cluster_cols=F, annotation_col=df)



