setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
library(plyr)

DB <- read.table("T1vsD7_DBpeaks_0812.info", header = F)
annotation <- read.csv("~/Desktop/Jiali/UTK/Apricot/Ppersica_298_v2.1.annotation_info.txt", header = T, stringsAsFactors = F, sep="\t",row.names = 1)
annotation <- annotation[!duplicated(annotation[,1]),] # remove the isoforms

head(DB)
names(DB)[5] <- "locusName"
DB$locusName <- gsub("\\.v2.1","",DB$locusName)
DB_annot <- join(DB,annotation[,c(1,10:12)],by="locusName", match="first")
write.csv(DB_annot,"T1vsD7_DBpeaks_0812.info.annot.csv")
