setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
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
