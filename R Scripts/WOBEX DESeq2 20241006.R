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
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, colData = metadata_filtered, design = ~ Pulldown + Run)
# Use varianceStabilizingTransformation for small datasets
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Check if the transformation worked correctly
plotPCA(vsd, intgroup = c("Run"))

# Perform variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Plot PCA, grouping by biological condition, pulldown, or run
plotPCA(vsd, intgroup = c("Pulldown", "Injury"))



