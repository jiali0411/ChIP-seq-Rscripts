setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
library(reshape2)
library(dplyr)
library(circlize)
#-----Circle plot----------------
#---------------------circal plot-----------------------------
## load chromosome info
cytoband.data <- read.table("~/Desktop/Jiali/UTK/Apricot/Pp_chrominfo.txt",header = F, sep = "\t") # peach chromosome range file, convert from gff
cytoband.df <- cytoband.data[,c(1,4,5)]
cytoband.df$V1 <- as.character(cytoband.df$V1)
cytoband.df$V4 <- as.numeric(cytoband.df$V4)
cytoband.df$V5 <- as.numeric(cytoband.df$V5)
par(mar=c(1,1,1,1))
## load ChIP peak data
peaks = read.csv("consensusBind.csv", header = T, row.names = 1)
peaks$CHR = paste0("Pp0",peaks$CHR)
# remove peaks not on the major 8 chromosomes 
keep <- c("Pp01", "Pp02", "Pp03", "Pp04", "Pp05", "Pp06", "Pp07", "Pp08")
filteredPeaks = peaks[peaks$CHR %in% keep, ]
T1 = filteredPeaks[,c(1,2,3,4)]
T1 <- filter(T1, X0 != 0)
T2 = filteredPeaks[,c(1,2,3,5)]
T2 <- filter(T2, X500 != 0)
T3 = filteredPeaks[,c(1,2,3,6)]
T3 <- filter(T3, X1000 != 0)
D3 = filteredPeaks[,c(1,2,3,7)]
D3 <- filter(D3, Day3 != 0)
D7 = filteredPeaks[,c(1,2,3,8)]
D7 <- filter(D7, Day7 != 0)
bed_list = list(T1,T2,T3,D3,D7)

# color map indicates the peak sites at each time point
circos.initializeWithIdeogram(cytoband.df, chromosome.index = "Pp01", plotType = c("axis", "labels"))
circos.genomicTrackPlotRegion(bed_list, stack = TRUE, panel.fun = function(region, value, ...) {
  i = getI(...)
  for(k in seq_len(nrow(region))) {
    circos.lines(rep(mean(region[k, 1], region[k, 2]), 2), c(i - 0.4, i + 0.4), straight = TRUE, col = "red")
  }
}, track.height = 0.2, bg.border = NA)
circos.clear()

# line graph
circos.initializeWithIdeogram(cytoband.df, plotType = c("axis", "labels"))
circos.genomicTrackPlotRegion(T1, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, area = TRUE, col = "blue", border = "blue", ...)
}, bg.border = NA, track.height = 0.08)
circos.genomicTrackPlotRegion(T2, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, area = TRUE, col = "blue", border = "blue", ...)
}, bg.border = NA, track.height = 0.08)
circos.genomicTrackPlotRegion(T3, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, area = TRUE, col = "blue", border = "blue", ...)
}, bg.border = NA, track.height = 0.08)
circos.genomicTrackPlotRegion(D3, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, area = TRUE, col = "blue", border = "blue", ...)
}, bg.border = NA, track.height = 0.08)
circos.genomicTrackPlotRegion(D7, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, area = TRUE, col = "blue", border = "blue", ...)
}, bg.border = NA, track.height = 0.08)
