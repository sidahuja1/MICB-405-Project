#Modified from From Lec 6-2

# Load libraries for data manipulation
library(tidyverse)
library(DESeq2)

# Read files into data frames
# Control replicate 1
sac_pombe_0h_rep1 <- read_tsv("0hr_R1ReadsPerGene.out.tab", # file name
                              col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                              skip = 4) # skip the first 4 lines
# Control replicate 2
sac_pombe_0h_rep2 <- read_tsv("0hr_R2ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)
# Control replicate 3
sac_pombe_0h_rep3 <- read_tsv("0hr_R3ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)
# Treatment replicate 1
sac_pombe_48h_rep1 <- read_tsv("48_rep1ReadsPerGene.out.tab",
                               col_names = c("gene_id", "total","antisense", "sense"),
                               skip = 4)

# Treatment replicate 2
sac_pombe_48h_rep2 <- read_tsv("48_rep2ReadsPerGene.out.tab",
                               col_names = c("gene_id", "total","antisense", "sense"),
                               skip = 4)
# Treatment replicate 3
sac_pombe_48h_rep3 <- read_tsv("48_rep3ReadsPerGene.out.tab",
                               col_names = c("gene_id", "total","antisense", "sense"),
                               skip = 4)

# DESeq2 is a widely used tool to assess truly differentially expressed genes
# Let's prepare the read counts for use in DESeq2

# Build a new data frame with gene names and read counts
# If you are using all stranded data, pick the appropriate column
# If using non-stranded data, use total column
sac_readcounts_0vs48 <- data.frame(row.names = sac_pombe_0h_rep1$gene_id,
                              sac_pombe_0h_rep1 = sac_pombe_0h_rep1$sense,
                              sac_pombe_0h_rep2 = sac_pombe_0h_rep2$sense,
                              sac_pombe_0h_rep3 = sac_pombe_0h_rep3$sense,
                              sac_pombe_48h_rep1 = sac_pombe_48h_rep1$sense,
                              sac_pombe_48h_rep2 = sac_pombe_48h_rep2$sense,
                              sac_pombe_48h_rep3 = sac_pombe_48h_rep3$sense)

# Clean the row names of an existing DESeq2 count data frame
rownames(sac_readcounts_0vs48) <- gsub("^.*_", "", rownames(sac_readcounts_0vs48))

# Let's check their order
colnames(sac_readcounts_0vs48)

# DESeq2 also requires the read counts to be in matrices, not data frames
# We can convert with as.matrix()
sac_matrix_0vs48 <- as.matrix(sac_readcounts_0vs48)

# DESeq2 also requires a table that provides information about the columns
# These must be in the order in which they appear on the count matrix
colnames(sac_readcounts_0vs48)

# Confusingly, column data is presented as a data frame The names of each sample
# should be the row name, and be identical and in order of appearance in the
# read counts matrix
columns_data_0vs48 <- data.frame(timepoint = c("0_hour",
                                         "0_hour",
                                         "0_hour",
                                         "48_hour",
                                         "48_hour",
                                         "48_hour"),
                           row.names = c("sac_pombe_0h_rep1",
                                         "sac_pombe_0h_rep2",
                                         "sac_pombe_0h_rep3",
                                         "sac_pombe_48h_rep1",
                                         "sac_pombe_48h_rep2",
                                         "sac_pombe_48h_rep3"))

# Create our DESeq2 object
dds_matrix_0vs48 <- DESeqDataSetFromMatrix(countData = sac_matrix_0vs48, #matrix 
                                     colData = columns_data_0vs48, #metadata file
                                     design = ~timepoint)
# Lets look at the object we created
dds_matrix_0vs48

# Set control condition using the relevel function
dds_matrix_0vs48$timepoint <- relevel(dds_matrix_0vs48$timepoint, ref = "0_hour")

# Run DESeq2
dds <- DESeq(dds_matrix_0vs48)

# These subsequent analyses are covered in Tutorial 7

# Log transforming count data
rld <- rlog(dds)
plotPCA(rld, intgroup = "timepoint")

# Distance matrix
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
colnames(sample_dist_matrix) <- NULL
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)

# names of the results that DESeq2 calculated
resultsNames(dds)

# Now we will extract the results for our comparison between the 12h timepoint
# and the 1h timepoint
res <- results(dds, name = "timepoint_48_hour_vs_0_hour") %>% as.data.frame()
# we save it as a data frame so we can apply dplyr functions better!
head(res)

# Getting rid of NA values
res_no_NA <- res %>% drop_na()
head(res_no_NA)
dim(res_no_NA) #look at how many rows you filtered out!

# Filtering for adjusted p-values <=0.05
res_filtered <- res_no_NA %>% filter(padj <= 0.05)
head(res_filtered)
dim(res_filtered) # look at how many rows you filtered out!

# Pull genes with more than 2x higher/lower expression
res_filtered_final <- res_filtered %>%
  filter(log2FoldChange <= -1 | log2FoldChange >= 1) %>%
rownames_to_column("gene_id")
# the '|' stand for OR here!
head(res_filtered_final)
dim(res_filtered_final)

# top 10 genes positively differentially expressed
top20_genes_0vs48 <- res_filtered_final %>%
  arrange(desc(log2FoldChange)) %>%
  # NOTE that we use the desc() function to organize the column in descending order
  head(n = 20)
top20_genes_0vs48

# how about the 10 genes negatively differentially expressed
bot20_genes_0vs48 <- res_filtered_final %>%
  arrange(log2FoldChange) %>%
  # NOTE since we don't use desc(), the column is organized in ascending order
  head(n = 20)
bot20_genes_0vs48

# Write CSV
write_csv(res_filtered_final, file = "results_filtered.csv")

# Separate filtered final into up and down reg genes
results_filtered_final_up <- res_filtered_final %>%
  filter(log2FoldChange >=1)

results_filtered_final_down <- res_filtered_final %>%
  filter(log2FoldChange <=-1)

#Write up and down reg genes into csv files for GO
write_csv(results_filtered_final_up, file = "results_filtered_up_0vs48.csv")
write_csv(results_filtered_final_down, file = "results_filtered_down_0vs48.csv")

# Load topGO
library(topGO)

# Read in output from DESeq2
up_gene <- read.csv("results_filtered_up_0vs48.csv")
down_gene <- read.csv("results_filtered_down_0vs48.csv")

# Read in a "gene universe" set
geneID2GO <- readMappings("gene_association_topGO.tsv")
geneUniverse <- names(geneID2GO)

# Get lists of upregulated/downregulated gene names
upregulated_gene_names <- as.character(up_gene$gene_id)
downregulated_gene_names <- as.character(down_gene$gene_id)

# Create a factor for upregulated genes
#up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_gene_names), 
#                       levels = c(0, 1),  # Ensure two levels: 0 = not upregulated, 1 = upregulated
#                       labels = c("Not Upregulated", "Upregulated"))
#names(up_gene_list) <- geneUniverse

# Create a factor for downregulated genes
#down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_gene_names),
#                         levels = c(0, 1),  # Ensure two levels: 0 = not downregulated, 1 = downregulated
#                         labels = c("Not Downregulated", "Downregulated"))
#names(down_gene_list) <- geneUniverse

#levels(up_gene_list)  # Should return: "Not Upregulated" and "Upregulated"
#levels(down_gene_list)  # Should return: "Not Downregulated" and "Downregulated"

# Factor the names
up_gene_list <- factor(as.integer(geneUniverse %in% upregulated_gene_names))
down_gene_list <- factor(as.integer(geneUniverse %in% downregulated_gene_names))
names(up_gene_list) <- geneUniverse
names(down_gene_list) <- geneUniverse

# Build a GOdata object in topGO for upregulated
up_GO_data <- new("topGOdata",
                  description = "sac_0h_48h",
                  ontology = "BP",
                  allGenes = up_gene_list,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)

down_GO_data <- new("topGOdata",
                    description = "sac_0h_48h",
                    ontology = "BP",
                    allGenes = down_gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# Run Fisher's test with the weight01 algorithm
up_result <- runTest(up_GO_data,
                     algorithm = "weight01",
                     statistic = "fisher")
down_result <- runTest(down_GO_data,
                       algorithm = "weight01",
                       statistic = "fisher")

# Summarize result for viewing
up_summary <- GenTable(up_GO_data,
                       weight01 = up_result,
                       orderBy = "up_result",
                       ranksOf = "up_result",
                       topNodes = 30)
down_summary <- GenTable(down_GO_data,
                         weight01 = down_result,
                         orderBy = "down_result",
                         ranksOf = "down_result",
                         topNodes = 30)

