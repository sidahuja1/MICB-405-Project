# Cen be modied for compariosn with 12 and 36 hour conditions

# Load libraries for data manipulation
suppressPackageStartupMessages(library(BiocManager))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
library(topGO)

# loading in the DESeq object
dds <- readRDS("dds.rds")

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

# Read in output from DESeq2
up_gene <- read.csv("results_filtered_up_0vs48.csv")
down_gene <- read.csv("results_filtered_down_0vs48.csv")

# Read in a "gene universe" set
geneID2GO <- readMappings("gene_association_topGO.tsv")
geneUniverse <- names(geneID2GO)

# Get lists of upregulated/downregulated gene names
upregulated_gene_names <- as.character(up_gene$gene_id)
downregulated_gene_names <- as.character(down_gene$gene_id)

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

