if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DiffBind")
#install.packages('httr')

library(DiffBind)
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

setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
Histone = dba(sampleSheet = "sample.csv")
dba.plotPCA(Histone, method = DBA_EDGER_BLOCK, attributes = c(DBA_TISSUE, DBA_REPLICATE))
# Rep1 samples cluster together and rep2 samples cluster together, maybe because of rep2 have no control
#------Find the consensus peaks found in both replicates-----------------
dba.plotVenn(Histone, Histone$masks$Day7)
Histone = dba.peakset(Histone, consensus = DBA_TISSUE, minOverlap = 2) # add consensus peaksets for peaks in both replicates
Histone_consensus = dba(Histone, mask = Histone$masks$Consensus)
Histone_consensus
dba.plotPCA(Histone_consensus, attributes = DBA_TISSUE, method = DBA_DESEQ2_BLOCK)


