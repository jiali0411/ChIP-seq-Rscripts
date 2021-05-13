setwd("~/Desktop/Jiali/UTK/Peach/USDA_peach/ChIP/ChIP-seq-Rscripts/")
library(rtracklayer)

## filter peaks not associated with a gene
annotation = import.gff3("../../Ppersica_298_v2.1.gene_exons.gff3") # gene positions
annotation_df <- data.frame(annotation)

gene_loc <- annotation_df[grep("gene",annotation_df$type),]
# add 1kb upstream and downstream
gene_loc$start <- gene_loc$start - 1000
gene_loc$end <- gene_loc$end + 1000
gene_range <- makeGRangesFromDataFrame(gene_loc, keep.extra.columns=TRUE,
                                    seqnames.field = "seqnames",
                                    start.field="start",
                                    end.field="end")
# overlap peaks and gene loci
overlapFilter <- function(filename){
peaks_all <- read.table(filename, header = F, sep = "\t")
peaks_range <- makeGRangesFromDataFrame(peaks_all, keep.extra.columns=TRUE,
                                      seqnames.field = "V1",
                                      start.field="V2",
                                      end.field="V3")
overlap_peaks <- subsetByOverlaps(peaks_range, gene_range, ignore.strand = TRUE)
filtered_peaks <- data.frame(overlap_peaks)
filtered_peaks <- filtered_peaks[,-c(4,5)]
write.table(filtered_peaks, paste0(filename,"_filtered"), quote = F, row.names = F, col.names = F, sep = "\t")
}

# run it on all files
file_list <- list.files(path = "./peaks",pattern = "*broadPeak", full.names = T)
for (file in file_list) {
  overlapFilter(file)
}
