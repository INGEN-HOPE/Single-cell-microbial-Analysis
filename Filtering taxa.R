#Load Library
library("phyloseq")
library("ggplot2")

#loading OTU table and taxa files
otu <- read.csv("OTU_table.csv", sep=",", row.names=1)
tax <- read.csv("taxa_table.csv", sep=",")
rownames(tax) <- rownames(otu)

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax))

#import sample information
sample <- read.csv("sample_table.csv", row.names=1)
sampledata = sample_data(sample)

#creating a phyloseq object
physeq1 = phyloseq(OTU, TAX, sampledata)

#agglomerating genus of same type
gp = tax_glom(physeq1, taxrank = "Genus")

#sorting in descending order
top <- names(sort(taxa_sums(gp), decreasing = TRUE))

#computing relative abundance
gp.prop <- transform_sample_counts(gp, function(x) x / sum(x))

#removing unwanted OTUs
gp.prop.top <- prune_taxa(top, gp.prop)

#combining OTU info, sample info and taxa annotations together
ps_df <- psmelt(gp.prop.top)

write.csv(ps_df,"Genus_level_abundance.csv")
