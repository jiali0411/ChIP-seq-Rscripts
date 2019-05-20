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
dba.plotPCA(Histone, method = DBA_EDGER_BLOCK, attributes = c(DBA_TISSUE, DBA_CALLER))
# Rep1 samples cluster together and rep2 samples cluster together, maybe because of rep2 have no control
#------Find the consensus peaks found in both replicates-----------------
# MACS2 and BCP have large variances, so analyse them seperately. 
dba.plotVenn(Histone, Histone$masks$`0` & Histone$masks$narrow)
Histone_macs = dba(Histone, mask = Histone$masks$narrow) # subset macs2 results
dba.plotHeatmap(Histone_macs) # heatmap of sample correlation 
Histone_macs = dba.peakset(Histone_macs, consensus=DBA_TISSUE, minOverlap = 2) # add consensus peaksets for peaks in both replicates

Histone_bcp = dba(Histone, mask = Histone$masks$bed) # subset BCP results
dba.plotHeatmap(Histone_bcp) # replicates have no correlation
Histone_bcp = dba.peakset(Histone_bcp, consensus = DBA_TISSUE, minOverlap = 2)
# retrive consensus peaks
Histone_macs_consensus = dba(Histone_macs, mask = Histone_macs$masks$Consensus) 
Histone_bcp_consensus = dba(Histone_bcp, mask=Histone_bcp$masks$Consensus)
dba.plotPCA(Histone_macs_consensus, attributes = c(DBA_TISSUE), method = DBA_EDGER)
dba.plotPCA(Histone_bcp_consensus, attributes = DBA_TISSUE, method = DBA_EDGER_BLOCK)

dba.plotHeatmap(Histone_macs_consensus, correlations=FALSE)
# Get binding scores for visualization
BindingScore = Histone_macs_consensus$binding
write.csv(BindingScore, "consensusBind.csv")

# Differential binding analysis
Histone_macs = dba.count(Histone_macs, minOverlap=2)
Histone_macs <- dba.contrast(Histone_macs, name1 = "chill", name2 = "warm", categories = DBA_TREATMENT)
Histone_macs <- dba.contrast(Histone_macs, name1 = "Endo", name2 = "PostEndo", group1 = Histone_macs$masks$Endo, group2 = Histone_macs$masks$postEndo)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "T2",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$`500`, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "T3",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$`1000`, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "D3",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$Day3, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T1", name2 = "D7",group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$Day7, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, name1 = 'T3', name2 = 'D3',group1 = Histone_macs$masks$`1000`, group2 = Histone_macs$masks$Day3)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T2", name2 = 'T3',group1 = Histone_macs$masks$`500`,group2 = Histone_macs$masks$`1000`)
Histone_macs <- dba.contrast(Histone_macs, name1 = "D3", name2 = "D7", group1 = Histone_macs$masks$Day3, group2 = Histone_macs$masks$Day7)
Histone_macs <- dba.contrast(Histone_macs, name1 = "T2", name2 = "D3", group1 = Histone_macs$masks$`500`, group2 = Histone_macs$masks$Day3)
# a contrast of endo vs postEndo automatically exists
dba.show(Histone_macs, bContrasts = T)

DffBind = dba.analyze(Histone_macs)
DffBind
DffBind_CvW = dba.report(DffBind, contrast = 1)
DffBind_EvP = dba.report(DffBind, contrast = 2)
DffBind_1v2 = dba.report(DffBind, contrast = 3)
DffBind_1v3 = dba.report(DffBind, contrast = 4)
DffBind_1vD3 = dba.report(DffBind, contrast = 5)
DffBind_1vD7 = dba.report(DffBind, contrast = 6)

corvals = dba.plotHeatmap(DffBind, contrast = 2, correlations = FALSE, attributes = DBA_TISSUE)

#-----------load gene info-----------------
genelocation <- read.table("~/Desktop/Jiali/UTK/Apricot/geneloc.txt", header = F)
genelocation <- genelocation[,c(1,4,5,9)]
genelocation$V9 <- colsplit(genelocation$V9,";",c('a', 'b'))
genelocation$V9$a <- gsub("ID=","",genelocation$V9$a)
genelocation$V9 <- genelocation$V9[,c(1)]
rownames(genelocation) <- genelocation$V9

write.table(DffBind_CvW, "chillvsWarm_DBpeaks.txt")
write.table(DffBind_EvP, "EndovsPost_DBpeaks.txt")
write.table(DffBind_1v2, "T1vsT2_DBpeaks.txt")
write.table(DffBind_1v3, "T1vsT3_DBpeaks.txt")
write.table(DffBind_1vD3, "T1vsD3_DBpeaks.txt")
write.table(DffBind_1vD7, "T1vsD7_DBpeaks.txt")
write.table(genelocation, "geneloc.txt")

