# Load library 
library(Seurat)

# Normalize the data in the "pathogen" assay
covid2 <- NormalizeData(covid, assay = "pathogen")

# Subset the Seurat object to include cells with more than 5 features in the "pathogen" assay
subset <- subset(covid2, subset = nFeature_pathogen > 5)

# Print the number of cells before subsetting
cat("Number of cells before subset:", ncol(covid2), "\n")

# Print the number of cells after subsetting
cat("Number of cells after subset:", ncol(subset), "\n")

# Print the number of cells removed during subsetting
cat("Number of cells removed:", ncol(covid2) - ncol(subset), "\n")

# Display the levels (cell identities)
levels(subset)
#
subset<- PrepSCTFindMarkers(subset)

# Find markers (differentially expressed genes) in the subsetted dataset
FindMarkers(subset, ident.1 = "Neutrophil_infected", ident.2 = "Neutrophil_recovered", logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1)

# Identify unique cell types in the subsetted Seurat object
cell_types <- unique(Idents(subset))

# Count the number of cells in each cell type
cell_type_counts <- table(Idents(subset))

# Print the cell type-wise count of cells
print(cell_type_counts)
