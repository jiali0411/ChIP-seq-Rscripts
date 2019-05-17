import sys

"""
Search the left closest/overlapped gene to the peak
Usage: python extractGenes.py <peak file> <gene locus file>
Output: overlapped gene as intragenic, non-overlapped gene as intergenic
"""

peakfile = sys.argv[1]
annotation = sys.argv[2]
# calculate closest gene to peak
def minDist(peak, annot):
    left = []
    for i in annot:
        if peak[0] == i[0]:
            dist = int(peak[1]) - int(i[2])
            left.append(dist)
            if dist <= 0:
                break
    closest = left.index(min(left))
    return annot[closest][3]

# read peak file
with open(peakfile, "r") as f:
    peaks = f.readlines()
peak_list =[]
for peak in peaks[1:]:
    position = peak.strip().split(" ")
    peak_list.append(position[1:4])
# read gene locus file
with open(annotation,"r") as a:
    loci = a.readlines()
# find closest gene at the peak
for j in peak_list:
    # print("is processing peak:",j)
    is_cover = False
    for i in loci[1:]:
        locus = i.strip().split(' ')
        is_cover = (j[0] == locus[1]) and (not (int(j[1]) > int(locus[3]) or int(j[2]) < int(locus[2])))
        if is_cover:
            print ' '.join(j),"intragenic",locus[4]  
            break 
    if (not is_cover):
        gene_list=[]
        for i in loci[1:]:
            locus = i.strip().split(' ')
            if j[0] == locus[1]:
                gene_list.append(locus[1:])
        closet_gene = minDist(j, gene_list)
        print ' '.join(j),"intergenic", closet_gene
    
