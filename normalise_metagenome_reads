
#Set working directory
getwd()
setwd(“<path/to/directory>”)

#Load library
library(metagenomeSeq)

#Upload the microbes count matrix  
raw_data<- read.csv(file ="count_matrix_microbes.csv", row.names = 1)

#create object
metaSeqObject <- newMRexperiment(raw_data) 

metaSeqObject_CSS <- cumNorm(metaSeqObject, p = cumNormStatFast(metaSeqObject))

OTU_read_count_CSS <- data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))

View(OTU_read_count_CSS)
write.csv(OTU_read_count_CSS,”count_matrix_normalized.csv")

