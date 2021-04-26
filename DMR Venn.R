setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")

# read diffbind files
DffBind_1v2 <- read.table("T1vsT2_DBpeaks_0812.txt", header = T)
DffBind_1v2$DMR <- paste0(DffBind_1v2$seqnames,"_",DffBind_1v2$start,":",DffBind_1v2$end)
DffBind_1v3 <- read.table("T1vsT3_DBpeaks_0812.txt", header = T)
DffBind_1v3$DMR <- paste0(DffBind_1v3$seqnames,"_",DffBind_1v3$start,":",DffBind_1v3$end)
DffBind_1vD3 <- read.table("T1vsD3_DBpeaks_0812.txt", header = T)
DffBind_1vD3$DMR <- paste0(DffBind_1vD3$seqnames,"_",DffBind_1vD3$start,":",DffBind_1vD3$end)
DffBind_1vD7 <- read.table("T1vsD7_DBpeaks_0812.txt", header = T)
DffBind_1vD7$DMR <- paste0(DffBind_1vD7$seqnames,"_",DffBind_1vD7$start,":",DffBind_1vD7$end)

# overlap DMR
peak_T1v2 <- DffBind_1v2$DMR
peak_T1v3 <- DffBind_1v3$DMR
peak_T1vD3 <- DffBind_1vD3$DMR
peak_T1vD7 <- DffBind_1vD7$DMR

# library
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Dark2")

#Make the plot
venn.diagram(
  x = list( peak_T1v2, peak_T1v3, peak_T1vD3, peak_T1vD7),
  category.names = c("T1vsT2" , "T1vsT3" , "T1vsD3","T1vsD7"),
  filename = 'plots/venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=myCol,
  fill = myCol,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer"
)

# heatmap
fold_change <- merge(DffBind_1v2[,c(9,12)], DffBind_1v3[,c(9,12)], by = "DMR", all=T)
colnames(fold_change) <- c("DMR","T2","T3")
fold_change <- merge(fold_change, DffBind_1vD3[,c(9,12)],by="DMR",all=T)
colnames(fold_change) <- c("DMR","T2","T3","D3")
fold_change <- merge(fold_change, DffBind_1vD7[,c(9,12)],by="DMR",all=T)
colnames(fold_change) <- c("DMR","T2","T3","D3","D7")
fold_change[is.na(fold_change)] <- 0
library(pheatmap)
pheatmap(-(fold_change[,-1]),scale = "none", cluster_cols = F,
         show_rownames = F, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         breaks = seq(-9, 9, length.out = 100))

