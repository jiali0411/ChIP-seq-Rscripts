setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
library(plyr)

DB <- read.table("T2vsT3_DBpeaks_2021.info", header = F)
annotation <- read.csv("~/Desktop/Jiali/UTK/Apricot/Ppersica_298_v2.1.annotation_info.txt", header = T, stringsAsFactors = F, sep="\t",row.names = 1)
annotation <- annotation[!duplicated(annotation[,1]),] # remove the isoforms

DE_data <- read.csv("../../analysis/DEG_T2vT3.csv", header = T)

head(DB)
names(DB)[13] <- "locusName"
names(DE_data)[1] <- "locusName"
# combine histone and DEG
DB_DE <- join(DB, DE_data[,c(1,3,7)], by="locusName", match="first")
# add gene annotation
DB_DE$locusName <- gsub("\\.v2.1","",DB_DE$locusName)
DB_annot <- join(DB_DE,annotation[,c(1,10:12)],by="locusName", match="first")
write.csv(DB_annot,"T2vsT3_DBpeaks_2021.info.annot.csv", row.names = F)
