# Load libraries
library(Seurat)
library(dplyr)
library(sctransform)
library(data.table)

#First, read the count matrix for human transcriptome, repeat for all samples
sample1_human_wta <- read.csv("path/to/dir/sample1_human_count_matrix.csv", row.names = 1)

#Next, read the count matrix for the AbSeq; repeat for all samples.
sample1_abseq_mat <- read.csv("path/to/dir/sample1_abseq_count_matrix.csv", row.names = 1)

#Finally, read the count matrix for the pathogens, repeat for all samples
sample1_microbe_mat <- read.csv("/path/to/dir/sample1_microbe_count_matrix.csv", row.names = 1)
#Create a SeuratObject using the human WTA counts, and add the AbSeq and pathogen counts as AssayObject. Repeat for all samples

sample1 <- CreateSeuratObject(counts = sample1_human_wta)
sample1_abseq <- CreateAssayObject(count = sample1_abseq_mat)
sample1_microbe <- CreateAssayObject(count = sample1_microbe_mat)
sample1 [["abseq"]] <- sample1_abseq
sample1 [["microbe"]] <- sample1_microbe

#Add labels defining the sample and group to each SeuratObject
sample1$sample <- "S1"
sample1$group <- "group1"

#Remove low quality cells, cells containing very low and high numbers of unique transcripts. Repeat for all samples
sample1 [["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^MT.")

#the cut-off parameters should be optimised based on the dataset 
sample1<- subset(sample1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20) 

#Merge all the sample level SeuratObject into one single object
covid <- merge(sample1, c(sample2,sample3,...), add.cell.ids = c("sample1", " sample2", " sample3",...))
Seurat::Project(object = covid) <- 'covid'

#Apply data normalization and batch effect correction
covid <- SplitObject(covid, split.by = " sample1")
covid <- lapply(X = covid, FUN = SCTransform)

#Perform the data integration
features <- SelectIntegrationFeatures(object.list = covid, nfeatures = 5000)
covid <- PrepSCTIntegration(object.list = covid, anchor.features = features)
covid.anchors <- FindIntegrationAnchors(object.list = covid, normalization.method = "SCT", anchor.features = features)
covid.sct <- IntegrateData(anchorset = covid.anchors, normalization.method = "SCT")
#Dimension Reduction
covid.sct <- RunPCA(covid.sct, verbose = TRUE)
DimPlot(covid, reduction="pca", raster = FALSE, group.by= "orig.ident") + NoLegend()
ElbowPlot(covid.sct) #Determine the dimensionality of the data
covid.int.sct <- RunTSNE(covid.int.sct, dims = 1:20, reduction = "pca")

#Identifying nearest neighbours and clustering of data
covid.int.sct <- FindNeighbors(covid.int.sct, reduction = "pca", dims = 1:20) 
covid.int.sct <- FindClusters(covid.int.sct, resolution = 0.4) 

#optimize the resolution according to the data
#Visualization of the clusters
DimPlot(covid.int.sct, reduction = "tsne", repel = TRUE, raster = FALSE)
DimPlot(covid.int.sct, reduction = "tsne", split.by = "sample", raster = FALSE)
DimPlot(covid.int.sct, reduction = "tsne", split.by = "group", raster = FALSE)

#Identify cluster-specific gene expression pattern
covidintmarkers <- FindAllMarkers(covid.int.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #min.pct and logfc.threshold to be optimized as per data

#Manually annotate the clusters based on cluster-specific gene expression pattern, and then add the cell type labels to the cluster
new.cluster.ids <- c("Naive T", "Naive B", "Classical Monocyte", "Neutrophil", "NK", "Memory B", "Macrophage", "Memory T", "Treg", "DC", "Platelet", "Plasma cell")

names(new.cluster.ids) <- levels(covid.int.sct)

covid.int.sct <- RenameIdents(covid.int.sct, new.cluster.ids)
DimPlot(covid.int.sct, reduction = "tsne", label = TRUE, pt.size = 0.5, raster = FALSE) + NoLegend()

#save the annotated count matrix as an RDS file
saveRDS(covid)

