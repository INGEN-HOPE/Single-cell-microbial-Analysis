# Load Libraries
library(data.table)
library(rtracklayer)

setwd("path/to/directory")

#extract all genes list from the genome gtf file
gtf_files <- list.files( pattern ="*.gtf")

for (gtf_file in gtf_files){
  # Import GTF file and convert to a data frame
  gtf.gr = rtracklayer::import(gtf_file)
  gtf.df = as.data.frame(gtf.gr)

  # Extract unique gene IDs
  genes = unique(gtf.df[ , "gene_id"])
  genes <-as.data.frame(genes)

 # Create output filename based on the GTF file
 output_file <- paste0(sub(".gtf$", "", basename(gtf_file)), "_genes.csv")

  write.csv(genes, file = output_file, row.names = FALSE)
print(paste("Processed:", gtf_file))
 }
