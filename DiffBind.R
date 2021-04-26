if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")
#install.packages('httr')

library(DiffBind)
#-------------------Test sample codes----------------------------
setwd(system.file("extra", package="DiffBind"))
tamoxifen = dba(sampleSheet = "tamoxifen.csv")
data("tamoxifen_counts")
data("tamoxifen_peaks")
tamoxifen = dba.count(tamoxifen, minOverlap = 3) # cannot access to read files
dba.overlap(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive, mode = DBA_OLAP_RATE)
dba.plotVenn(tamoxifen, tamoxifen$masks$MCF7 & tamoxifen$masks$Responsive)

tamoxifen = dba.peakset(tamoxifen, consensus = c(DBA_TISSUE, DBA_CONDITION), minOverlap = 0.66)
tamoxifen_consensus = dba(tamoxifen, mask = tamoxifen$masks$Consensus)
tamoxifen_consensus
tamoxifen = dba.count(tamoxifen, peaks = tamoxifen$masks$Consensus) # won't work because no bam files are available in this package.
dba.plotPCA(tamoxifen, th=0.05)
tamoxifen.DB = dba.report(tamoxifen)
corvals = dba.plotHeatmap(tamoxifen, correlations = FALSE)
dba.plotPCA(tamoxifen, method = DBA_EDGER_BLOCK, attributes = c(DBA_TISSUE, DBA_CONDITION, DBA_REPLICATE))


sample <- read.csv("tamoxifen.csv", header = T) 
# I follow the tamoxifen.csv file format to create my own data, in the working folder as sample.csv
#-----------------load my own data-------------------------------
setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
Histone = dba(sampleSheet = "sample.csv")
pcaPlot <- dba.plotPCA(Histone, method = DBA_EDGER_BLOCK, attributes = c(DBA_TISSUE))
# Rep1 samples cluster together and rep2 samples cluster together, maybe because of rep2 have no control
#------Find the consensus peaks found in both replicates-----------------
dba.plotHeatmap(Histone, ColAttributes = DBA_TREATMENT) # heatmap of sample correlation 
Histone_macs = dba.peakset(Histone, consensus=DBA_TISSUE, minOverlap = 2) # add consensus peaksets for peaks in both replicates
Histone_macs

# correlation scores
Histone.corRec <- dba.overlap(Histone_macs,mode=DBA_OLAP_ALL,bCorOnly=FALSE)

# retrive consensus peaks
Histone_macs_consensus = dba(Histone_macs, mask = Histone_macs$masks$Consensus)
dba.plotPCA(Histone_macs_consensus, attributes = c(DBA_TISSUE), method = DBA_EDGER)
dba.plotHeatmap(Histone_macs_consensus, ColAttributes = DBA_TREATMENT,correlations=FALSE)
# Get binding scores for visualization
BindingScore = Histone_macs_consensus$binding
write.csv(BindingScore, "consensusBind.csv")
# separate bindscore for each timepoint
T1 <- data.frame(Histone_macs_consensus$peaks[[1]])
T1 <- T1[grep("Pp",T1$V1),] # filter out small scaffolds
T1$V4 <- T1$V4*100
write.table(T1, "consensusPeaks_T1.txt", col.names = F, row.names = F, sep = "\t", quote = F)

T2 <- data.frame(Histone_macs_consensus$peaks[[2]])
T2 <- T2[grep("Pp",T2$V1),] # filter out small scaffolds
T2$V4 <- T2$V4 * 100
write.table(T2, "consensusPeaks_T2.txt", col.names = F, row.names = F, sep = "\t", quote = F)

T3 <- data.frame(Histone_macs_consensus$peaks[[3]])
T3 <- T3[grep("Pp",T3$V1),] # filter out small scaffolds
T3$V4 <- T3$V4*100
write.table(T3, "consensusPeaks_T3.txt", col.names = F, row.names = F, sep = "\t", quote = F)

D3 <- data.frame(Histone_macs_consensus$peaks[[4]])
D3 <- D3[grep("Pp",D3$V1),] # filter out small scaffolds
D3$V4 <- D3$V4*100
write.table(D3, "consensusPeaks_D3.txt", col.names = F, row.names = F, sep = "\t", quote = F)

D7 <- data.frame(Histone_macs_consensus$peaks[[5]])
D7 <- D7[grep("Pp",D7$V1),] # filter out small scaffolds
D7$V4 <- D7$V4*100
write.table(D7, "consensusPeaks_D7.txt", col.names = F, row.names = F, sep = "\t", quote = F)

# Differential binding analysis
Histone_macs = dba.count(Histone, minOverlap=2)
Histone_macs <- dba.contrast(Histone_macs, name1 = "chill", name2 = "warm", categories = DBA_TREATMENT)
Histone_macs <- dba.contrast(Histone_macs, name1 = "Endo", name2 = "PostEndo", group1 = Histone_macs$masks$Endo, group2 = Histone_macs$masks$postEndo)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "T2",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$`500`, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "T3",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$`1000`, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "D3",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$Day3, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "D7",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$Day7, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = 'T3', name2 = 'D3',group1 = Histone_macs$masks$`1000`, group2 = Histone_macs$masks$Day3)
Histone_macs <- dba.contrast(Histone_macs, name1 = 'T3', name2 = 'D7',group1 = Histone_macs$masks$`1000`, group2 = Histone_macs$masks$Day7)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T2", name2 = 'T3',group1 = Histone_macs$masks$`500`,group2 = Histone_macs$masks$`1000`)
Histone_macs <- dba.contrast(Histone_macs, name1 = "D3", name2 = "D7", group1 = Histone_macs$masks$Day3, group2 = Histone_macs$masks$Day7)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T2", name2 = "D3", group1 = Histone_macs$masks$`500`, group2 = Histone_macs$masks$Day3)
# a contrast of endo vs postEndo automatically exists
dba.show(Histone_macs, bContrasts = T)

DffBind = dba.analyze(Histone_macs)
DffBind
DffBind_CvW = dba.report(DffBind, contrast = 1)
DffBind_EvP = dba.report(DffBind, contrast = 2)
DffBind_1v2 = dba.report(DffBind, contrast = 1)
DffBind_1v3 = dba.report(DffBind, contrast = 2)
DffBind_1vD3 = dba.report(DffBind, contrast = 3)
DffBind_1vD7 = dba.report(DffBind, contrast = 4)

#-----------load gene info-----------------
library(reshape2)
genelocation <- read.table("~/Desktop/Jiali/UTK/Apricot/geneloc.txt", header = F)
genelocation <- genelocation[,c(1,4,5,9)]
genelocation$V9 <- colsplit(genelocation$V9,";",c('a', 'b'))
genelocation$V9$a <- gsub("ID=","",genelocation$V9$a)
genelocation$V9 <- genelocation$V9[,c(1)]
rownames(genelocation) <- genelocation$V9

write.table(DffBind_CvW, "chillvsWarm_DBpeaks.txt")
write.table(DffBind_EvP, "EndovsPost_DBpeaks.txt")
write.table(DffBind_1v2, "T1vsT2_DBpeaks_0812.txt")
write.table(DffBind_1v3, "T1vsT3_DBpeaks_0812.txt")
write.table(DffBind_1vD3, "T1vsD3_DBpeaks_0812.txt")
write.table(DffBind_1vD7, "T1vsD7_DBpeaks_0812.txt")
write.table(genelocation, "geneloc.txt")

# ------------compare T3 with D3 and D7--------------------
DffBind_3vD3 = dba.report(DffBind, contrast = 1)
DffBind_3vD7 = dba.report(DffBind, contrast = 2)
write.table(DffBind_3vD7, "T3vsD7_DBpeaks.txt")

