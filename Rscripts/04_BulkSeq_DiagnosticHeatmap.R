rm(list=ls())

##Here, we preform determinstic clustering on sample-to-sample distances

#First, we load our Deseq2 object as it contains all the counts and metadata
load('R_objects/dds_BulkSeq_Aging.bin')

dds_AllRegions <- dds_CA1_list$All

#We'll use Deseq2's vst function to perform varianceStabilizingTransformation
vsd <- vst(dds_AllRegions)
#We'll focus on the 20,000 highest expressed genes, so we'll select those
select <- order(rowMeans(counts(dds_AllRegions,normalized=F)),decreasing=TRUE)[1:20000]
#Calculate sample-to-sample distance matrix and plot heatmap upon h-clustering
distsVSD <- dist(t(assay(vsd)[select,]))     
#Let's reformat that as a matrix
mat <- as.matrix(distsVSD)

#We want use the ward.D2 method for clustering so we'll have to define that using hclust
hc <- hclust(distsVSD, method = 'ward.D2')

#####Now we plot Figure S2E
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none", col = colorRampPalette(c("red",'yellow', "blue")) , margin=c(13, 13)
