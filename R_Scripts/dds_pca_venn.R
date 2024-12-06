library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)

# Loading in STAR output files into R
sp_0h_rep1 <- read_tsv("data/00hr_R1ReadsPerGene.out.tab", # file name
                              col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                              skip = 4) # skip the first 4 lines

sp_0h_rep2 <- read_tsv("data/00hr_R2ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_0h_rep3 <- read_tsv("data/00hr_R3ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_12h_rep1 <- read_tsv("data/12hr_R1ReadsPerGene.out.tab", # file name
                              col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                              skip = 4) # skip the first 4 lines

sp_12h_rep2 <- read_tsv("data/12hr_R2ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_12h_rep3 <- read_tsv("data/12hr_R3ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_24h_rep1 <- read_tsv("data/24hr_R1ReadsPerGene.out.tab", # file name
                              col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                              skip = 4) # skip the first 4 lines

sp_24h_rep2 <- read_tsv("data/24hr_R2ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_24h_rep3 <- read_tsv("data/24hr_R3ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_36h_rep1 <- read_tsv("data/36hr_R1ReadsPerGene.out.tab", # file name
                              col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                              skip = 4) # skip the first 4 lines

sp_36h_rep2 <- read_tsv("data/36hr_R2ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_36h_rep3 <- read_tsv("data/36hr_R3ReadsPerGene.out.tab",
                              col_names = c("gene_id", "total","antisense", "sense"),
                              skip = 4)

sp_48h_rep1 <- read_tsv("data/48hr_R1ReadsPerGene.out.tab", # file name
                        col_names = c("gene_id", "total", "antisense", "sense"), # name the columns
                        skip = 4) # skip the first 4 lines

sp_48h_rep2 <- read_tsv("data/48hr_R2ReadsPerGene.out.tab",
                        col_names = c("gene_id", "total","antisense", "sense"),
                        skip = 4)

sp_48h_rep3 <- read_tsv("data/48hr_R3ReadsPerGene.out.tab",
                        col_names = c("gene_id", "total","antisense", "sense"),
                        skip = 4)

# Creating the counts matrix
dat <- data.frame(row.names = sp_0h_rep1$gene_id,
                  sp_0h_rep1 = sp_0h_rep1$sense,
                  sp_0h_rep2 = sp_0h_rep2$sense,
                  sp_0h_rep3 = sp_0h_rep3$sense,
                  sp_12h_rep1 = sp_12h_rep1$sense,
                  sp_12h_rep2 = sp_12h_rep2$sense,
                  sp_12h_rep3 = sp_12h_rep3$sense,
                  sp_24h_rep1 = sp_24h_rep1$sense,
                  sp_24h_rep2 = sp_24h_rep2$sense,
                  sp_24h_rep3 = sp_24h_rep3$sense,
                  sp_36h_rep1 = sp_36h_rep1$sense,
                  sp_36h_rep2 = sp_36h_rep2$sense,
                  sp_36h_rep3 = sp_36h_rep3$sense,
                  sp_48h_rep1 = sp_48h_rep1$sense,
                  sp_48h_rep2 = sp_48h_rep2$sense,
                  sp_48h_rep3 = sp_48h_rep3$sense)

dat_matrix<- as.matrix(dat) 

# Creating the metadata data frame
metadata <- data.frame(row.names = colnames(dat_matrix), 
                       Timepoint = c("0_hour", "0_hour", "0_hour",
                                     "12_hour", "12_hour", "12_hour",
                                     "24_hour", "24_hour", "24_hour",
                                     "36_hour", "36_hour", "36_hour",
                                     "48_hour", "48_hour", "48_hour"))

# Creating DESeq object and setting 0_hour as control
dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix,
                                     colData = metadata,
                                     design = ~Timepoint)

dds_matrix$Timepoint <- relevel(dds_matrix$Timepoint, ref = "0_hour")

dds <- DESeq(dds_matrix)

saveRDS(dds, "dds.rds")

