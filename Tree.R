# Caitlin Falconer 22/10/17
# Plot and save dendrogram based on distance matrix
DIR = commandArgs(TRUE)
INFILE=paste(DIR, 'distmat.csv', sep='/')
OUTFILE=paste(DIR, 'LargeSVtree.jpeg', sep='/')
data = read.csv(INFILE, header=T, sep="\t")
dmat = dist(data, method='binary')
clust = hclust(dmat)
jpeg(file=OUTFILE)
  plot(clust, main="Large Structural Variant Phylogenetic Tree", xlab="dist(distmat.csv, method='binary')")
dev.off()
