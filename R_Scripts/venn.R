library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(VennDiagram)

# loading in the DESeq object
dds <- readRDS("dds.rds")

# Filtering out ncRNA genes
ncrna_genes <- grepl("NCRNA", rownames(dds), ignore.case = F)

dds_filtered <- dds[!ncrna_genes, ]

# Removing SPOM_ prefix from gene IDs
rownames(dds_filtered) <- gsub("SPOM_", "", rownames(dds_filtered))

# Pairwise comparisons (DGE)
res_0_vs_12 <- results(dds_filtered, contrast=c("Timepoint", "0_hour", "12_hour"))
res_0_vs_24 <- results(dds_filtered, contrast=c("Timepoint", "0_hour", "24_hour"))
res_0_vs_36 <- results(dds_filtered, contrast=c("Timepoint", "0_hour", "36_hour"))
res_0_vs_48 <- results(dds_filtered, contrast=c("Timepoint", "0_hour", "48_hour"))

# Function to filter over expressed genes in the 0 hour vs other timepoints
filter_genes <- function(res) {
  res_filtered <- res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >= 2, ]
  return(rownames(res_filtered))
}

# Running function
overexpressed_0_vs_12 <- filter_genes(res_0_vs_12)
overexpressed_0_vs_24 <- filter_genes(res_0_vs_24)
overexpressed_0_vs_36 <- filter_genes(res_0_vs_36)
overexpressed_0_vs_48 <- filter_genes(res_0_vs_48)

# Combining the genes
gene_lists <- list(
  `0_vs_12` = overexpressed_0_vs_12,
  `0_vs_24` = overexpressed_0_vs_24,
  `0_vs_36` = overexpressed_0_vs_36,
  `0_vs_48` = overexpressed_0_vs_48
)

# Selecting color pallete
myCol <- brewer.pal(4, "Set1")

# Creating venn diagram
venn.plot <- draw.quad.venn(
  area1 = length(overexpressed_0_vs_12),
  area2 = length(overexpressed_0_vs_24),
  area3 = length(overexpressed_0_vs_36),
  area4 = length(overexpressed_0_vs_48),
  n12 = length(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24)),
  n13 = length(intersect(overexpressed_0_vs_12, overexpressed_0_vs_36)),
  n14 = length(intersect(overexpressed_0_vs_12, overexpressed_0_vs_48)),
  n23 = length(intersect(overexpressed_0_vs_24, overexpressed_0_vs_36)),
  n24 = length(intersect(overexpressed_0_vs_24, overexpressed_0_vs_48)),
  n34 = length(intersect(overexpressed_0_vs_36, overexpressed_0_vs_48)),
  n123 = length(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24), overexpressed_0_vs_36)),
  n124 = length(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24), overexpressed_0_vs_48)),
  n134 = length(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_36), overexpressed_0_vs_48)),
  n234 = length(intersect(intersect(overexpressed_0_vs_24, overexpressed_0_vs_36), overexpressed_0_vs_48)),
  n1234 = length(intersect(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24), overexpressed_0_vs_36), overexpressed_0_vs_48)),
  category = c("12 hour", "24 hour", "36 hour", "48 hour"),
  fill = myCol,
  lty = "blank",
  fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cex = 1.25,
  cat.cex = 1.4
)


# Function to filter under expressed genes in the 0 hour vs other timepoints
filter_genes <- function(res) {
  res_filtered <- res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > -2, ]
  return(rownames(res_filtered))
}

# Running function
overexpressed_0_vs_12 <- filter_genes(res_0_vs_12)
overexpressed_0_vs_24 <- filter_genes(res_0_vs_24)
overexpressed_0_vs_36 <- filter_genes(res_0_vs_36)
overexpressed_0_vs_48 <- filter_genes(res_0_vs_48)

# Combining the genes
gene_lists <- list(
  `0_vs_12` = overexpressed_0_vs_12,
  `0_vs_24` = overexpressed_0_vs_24,
  `0_vs_36` = overexpressed_0_vs_36,
  `0_vs_48` = overexpressed_0_vs_48
)

# Creating venn diagram
venn.plot <- draw.quad.venn(
  area1 = length(overexpressed_0_vs_12),
  area2 = length(overexpressed_0_vs_24),
  area3 = length(overexpressed_0_vs_36),
  area4 = length(overexpressed_0_vs_48),
  n12 = length(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24)),
  n13 = length(intersect(overexpressed_0_vs_12, overexpressed_0_vs_36)),
  n14 = length(intersect(overexpressed_0_vs_12, overexpressed_0_vs_48)),
  n23 = length(intersect(overexpressed_0_vs_24, overexpressed_0_vs_36)),
  n24 = length(intersect(overexpressed_0_vs_24, overexpressed_0_vs_48)),
  n34 = length(intersect(overexpressed_0_vs_36, overexpressed_0_vs_48)),
  n123 = length(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24), overexpressed_0_vs_36)),
  n124 = length(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24), overexpressed_0_vs_48)),
  n134 = length(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_36), overexpressed_0_vs_48)),
  n234 = length(intersect(intersect(overexpressed_0_vs_24, overexpressed_0_vs_36), overexpressed_0_vs_48)),
  n1234 = length(intersect(intersect(intersect(overexpressed_0_vs_12, overexpressed_0_vs_24), overexpressed_0_vs_36), overexpressed_0_vs_48)),
  category = c("12 hour", "24 hour", "36 hour", "48 hour"),
  fill = myCol,
  lty = "blank",
  fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cex = 1.25,
  cat.cex = 1.4
)
