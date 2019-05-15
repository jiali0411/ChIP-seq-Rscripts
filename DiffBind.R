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
#-----------------load own data-------------------------------
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
Histone_macs <- dba.contrast(Histone_macs, group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$`500`, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$`1000`, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$Day3, categories=DBA_TISSUE)
Histone_macs <- dba.contrast(Histone_macs, group1 = Histone_macs$masks$`0`, group2 = Histone_macs$masks$Day7, categories=DBA_TISSUE)
dba.show(Histone_macs, bContrasts = T)

DffBind = dba.analyze(Histone_macs)

