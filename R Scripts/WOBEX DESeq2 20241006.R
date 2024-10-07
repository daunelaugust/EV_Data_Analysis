# Load necessary libraries
library(DESeq2)
library(tidyverse)


# save filenames in vars
meta <- paste(getwd(),"/Excel Data Files/Master WOBEX Spreadsheet 20241006.csv",sep ='')
counts <- (paste(getwd(),"/Excel Data Files/WOBEX ReadCountMatrix 20241006.csv",sep =''))

# Reading metadata and count matrix
metadata <- read.csv(meta,row.names =1 )
count_matrix <- as.matrix(read.csv(counts, fileEncoding="UTF-8-BOM",row.names=1))


head(metadata)
head(count_matrix)

# making sure the row names in metadata matches to column names in count_matrix
all(colnames(count_matrix) %in% rownames(metadata))

# are they in the same order?
all(colnames(count_matrix) == rownames(metadata))


# filter Samples with Less Than 100 Reads
sample_sums <- colSums(count_matrix)
filtered_samples <- sample_sums[sample_sums >= 100]  # keeps samples with >= 100 reads
count_matrix_filtered <- count_matrix[, names(filtered_samples)]

# filter miRNAs with Less Than 10 Reads
miRNA_sums <- rowSums(count_matrix_filtered)
filtered_miRNAs <- miRNA_sums[miRNA_sums >= 10]  # keeps miRNAs with >= 10 reads
count_matrix_filtered <- count_matrix_filtered[rownames(count_matrix_filtered) %in% names(filtered_miRNAs), ]

# Filter metadata to match the remaining samples
metadata_filtered <- metadata[rownames(metadata) %in% colnames(count_matrix_filtered), ]

# Ensure metadata is in the same order as count_matrix
metadata_filtered <- metadata_filtered[colnames(count_matrix_filtered), ]

# Another check to confirm that the column names of counts match row names of metadata
all(colnames(count_matrix_filtered) %in% rownames(metadata_filtered)) 
all(colnames(count_matrix_filtered) == rownames(metadata_filtered))  



# Transform the count data for PCA
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, colData = metadata_filtered, design = ~Run * Pulldown)

# Use varianceStabilizingTransformation for small datasets
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Check if the transformation worked correctly
plotPCA(vsd, intgroup = c("Run"))


dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, colData = metadata_filtered, design = ~ Run + Pulldown + Injury + Timepoint)
dds <- DESeq(dds)


# Model with interactions between Injury, Timepoint, and Pulldown
dds2 <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, 
                              colData = metadata_filtered, 
                              design = ~ Run + Pulldown * Injury * Timepoint)

# Running the DESeq2 analysis
dds2 <- DESeq(dds2)

# Subset for GLAST and pre-injury samples
dds_glast_pre_injury <- dds[ , dds$Timepoint == "PRE" & dds$Pulldown == "GLAST"]

# Extract results for GLAST pre-injury (Injured vs Sham)
res_glast_pre_injury <- results(dds_glast_pre_injury, contrast = c("Injury", "Injured", "Sham"))

# Summary of results

# Subset for GluR2 and pre-injury samples
dds_glur2_pre_injury <- dds[ , dds$Timepoint == "PRE" & dds$Pulldown == "GluR2"]

# Extract results for GluR2 pre-injury (Injured vs Sham)
res_glur2_pre_injury <- results(dds_glur2_pre_injury, contrast = c("Injury", "Injured", "Sham"))

# Summary of results
summary(res_glur2_pre_injury)

summary(res_glast_pre_injury)

# Filter for significant miRNAs (padj < 0.05)
sig_mirna_glast <- res_glast_pre_injury[which(res_glast_pre_injury$padj < 0.05), ]

# Show the names of the significant miRNAs
sig_mirna_names_glast <- rownames(sig_mirna_glast)
print(sig_mirna_names_glast)

sig_mirna_glur2 <- res_glur2_pre_injury[which(res_glur2_pre_injury$padj < 0.05), ]

# Show the names of the significant miRNAs for GluR2
sig_mirna_names_glur2 <- rownames(sig_mirna_glur2)
print(sig_mirna_names_glur2)



# Apply the variance-stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extract PCA data from vsd
pcaData <- plotPCA(vsd, intgroup = "Run", returnData = TRUE)

# Plot using ggplot2
library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color = Run)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot of Samples Colored by Sequencing Run") +
  xlab(paste0("PC1: ", round(100 * attr(pcaData, "percentVar")[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pcaData, "percentVar")[2], 1), "% variance"))


