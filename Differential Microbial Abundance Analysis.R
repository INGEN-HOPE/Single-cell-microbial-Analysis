# Load libraries
library(tidyverse)
library(cowplot)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(Matrix.utils) 

#loading count data
counts<-read.csv("microbial_counts_samplewise.csv", row.names = "X")
counts<-counts[,c(-1,-2)]
t_counts<-t(counts)
 
#loading metadata
metadata<-read.csv("microbial_counts_metadata.csv", row.names = "X")
metadata[,'celltype']<-as.factor(metadata[,'celltype'])
metadata[,'sample']<-as.factor(metadata[,'sample'])
metadata[,'group']<-as.factor(metadata[,'group'])
str(metadata)

#Create SingleCellExperiment object
sce<-SingleCellExperiment(assays= list(counts =counts), colData=metadata)
#Extract unique cluster and sample names
cluster_names <- levels(colData(sce)$celltype)
sample_names <- levels(colData(sce)$sample)
 
#Subset metadata for aggregation
groups <- colData(sce)[, c("celltype", "sample")]
head(groups)

#Aggregating counts across cluster-sample groups
#transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)),  groupings = groups, fun = "sum")
aggr_counts <- t(aggr_counts)
 
#Exploring structure of function output
tstrsplit(colnames(aggr_counts), "_") %>% str()
 
#Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)

#Extract indices of columns matching a specific cell type
b_cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "Naive B")
b_cell_idx
 
#Loop over all cell types to extract corresponding counts, and store information in a list
#Initialise list for counts by cluster
counts_ls <- list()
for (i in seq_along(cluster_names)) {
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]]      == cluster_names[i])
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}

#Explore the different components of the list
str(counts_ls)
 
#Prepare metadata for DE analysis
metadata <- colData(sce) %>%
          as.data.frame() %>%
          dplyr::select(group, sample)
 
 
#Exclude duplicate rows
metadata <- metadata[!duplicated(metadata), ]
 
rownames(metadata) <- metadata$sample
head(metadata)

# Number of cells per sample and cluster
t <- table(colData(sce)$sample, colData(sce)$celltype)
 
#Initialise list for metadata by cluster
metadata_ls <- list()
 
for (i in 1:length(counts_ls)) {
 
#Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
 
#Use tstrsplit() to separate cluster (cell type) and sample IDs
df$celltype <- tstrsplit(df$cluster_sample_id, "_")[[1]]
df$batch <- tstrsplit(df$cluster_sample_id, "_")[[2]]
df$sample  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
df$sample_id<-paste0(df$batch,"_",df$sample)
 
df$group<-metadata$group[match((df$sample_id), row.names(metadata))]

  #Retrieve cell count information for this cluster from global cell count table
idx <- which(colnames(t) == unique(df$celltype))
cell_counts <- t[, idx]
 
#Remove samples with zero cell contributing to the cluster
cell_counts <- cell_counts[cell_counts > 0]
 
#Match order of cell_counts and sample_ids
sample_order <- match(df$sample_id, names(cell_counts))
cell_counts <- cell_counts[sample_order]
 
#Append cell_counts to data frame
df$cell_count <- cell_counts

#Join data frame (capturing metadata specific to cluster) to generic metadata
df <- plyr::join(df, metadata,
               	by = intersect(names(df), names(metadata)))
 
#Update rownames of metadata to match colnames of count matrix, as needed later for DE
rownames(df) <- df$cluster_sample_id
 
#Store complete metadata for cluster i in list
metadata_ls[[i]] <- df
names(metadata_ls)[i] <- unique(df$celltype)
}
 
#Explore the different components of the list
str(metadata_ls)

#Starting pseudo bulk analysis using DESEq2#
 
# Select cell type of interest
cluster_names
 
# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))
 
idx <- which(names(counts_ls) == "DC")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]
 
# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))
 
cluster_counts<-as.matrix(cluster_counts +1)#adding 1 to the count matrix to allow log transformation of the data

# Create DESeq2 object       
dds <- DESeqDataSetFromMatrix(cluster_counts,
                              colData = cluster_metadata,
                              design = ~ group)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
 
# Plot PCA
DESeq2::plotPCA(rld, intgroup = "cell_count")
 
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
 
# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
res <-results(dds)




