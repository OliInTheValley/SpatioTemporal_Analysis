rm(list=ls())

##Here, we perform an analysis in Seurat, akin to single-cell data to identify the 'purity' of our bulk tissue RNA-seq data
#Then we'll perform marker analysis to develop signatures for each region that we can plot onto spatial transcriptome data

#First, we load our Deseq2 object as it contains all the counts and metadata
load('R_objects/dds_BulkSeq_Aging.bin')

#We create a Seurat object using the counts from the 'All' Desea2 object list (that's where we stored data from all the tissues)
CA1_seurat <- CreateSeuratObject(counts(dds_CA1_list$All))
#We'll also carry over some of the metadata
CA1_seurat@meta.data[,4:9] <- as.data.frame(colData(dds_CA1_list$All))[c('tissueLong', 'sex', 'age', 'mouseID', 'tissue')]
#We perform standard Log-normalization
CA1_seurat <- NormalizeData(CA1_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#We define the variable genes (2000 should do)
CA1_seurat <- FindVariableFeatures(CA1_seurat, selection.method = "vst", nfeatures = 2000)

#We'll inspect the most variable features manually
#first we'll define the top 10 variable features
top10 <- head(VariableFeatures(CA1_seurat), 10)
#And then use Seurat variablefeatureplot function
plot1 <- VariableFeaturePlot(CA1_seurat)
#Let's label these points using the labelpoints function
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Now we'll scale the data using all genes in the count matrix as input
CA1_seurat <- ScaleData(CA1_seurat, features = rownames(CA1_seurat), verbose = T)

#Now we'll run the PCA and visualize the explained variance per PC
CA1_seurat <- RunPCA(CA1_seurat)
ElbowPlot(CA1_seurat, ndims = 40)

#Now we'll run the graph-based clustering across the first 40 PCs
CA1_seurat <- FindNeighbors(CA1_seurat, dims = 1:40)
CA1_seurat <- FindClusters(CA1_seurat, resolution = 0.8)
#For visualization, we'll also run a UMAP across the same parameters
CA1_seurat <- RunUMAP(CA1_seurat, dims = 1:40, verbose = T)
#We'll set the default idents to be the same as in the 'tissue' column
Idents(CA1_seurat) <- 'tissue'

#Let's visualise the outcome 
#####And plot Figure 1C
DimPlot(CA1_seurat, group.by = 'tissue', label = T)
#####plot Figure S1B
DimPlot(CA1_seurat, group.by = 'sex', cols = c('blue', 'red'), label = T)
#####plot Figure S1C
DimPlot(CA1_seurat, group.by = 'age', label = T)
#####plot Figure S1D
DimPlot(CA1_seurat, group.by = 'seurat_clusters', label = T)

#We'll save the object for alter usage
saveRDS(CA1_seurat, file = 'R_objects/CA1_seurat.rds')

###Now we'll run a marker analysis which we'll use to build signatures of region-enriched genes
#Since this is bulk data, we will use the Deseq2 method that's built into seurat
tissue_markers <- FindAllMarkers(CA1_seurat, only.pos = T, verbose = T, test.use = 'DESeq2')

#We'll write out the marker lists
write.table(tissue_markers, quote = F, sep = '\t', col.names = NA, file='Output_tables/CA1_brainregionMarkers_Deseq2.txt')

#Now we'll use VISION's workflow to create signatures based on all signinficant marker genes
#We'll loop over every region/tissue in our marker list, extract the significant marker genes and build a signature vector
#First, we'll set up the signature vector
mySignatures_RegionMarkersDeseq2 <- c()
for (tissue in unique(tissue_markers$cluster)) {
  #Extract the respetive region's significant marker genes
  marker_genes <- na.omit((tissue_markers %>% filter(cluster==tissue) %>% filter(p_val_adj < 0.05))$gene)
  #VISION requires a named vector where the values are -1 or 1 (for up/down regulation); so we'll create a vector with +1 of the same length as the marker genes
  signature_vector <- rep(1, length(marker_genes))
  #Now we'll name that vector. Vision requirs all genes to be uppercase
  names(signature_vector) <- toupper(marker_genes)
  #We'll remove any potential dupplicates that could cause trouble down the line
  signature_vector <- signature_vector[!(duplicated(names(signature_vector)) | duplicated(names(signature_vector),fromLast=TRUE))]
  
  sig_tissue_markergenes <- createGeneSignature(name = paste("CA1_", tissue, "_RegionMarkersDeseq2", sep = ''), sigData = signature_vector)
  mySignatures_RegionMarkersDeseq2 <- c(mySignatures_RegionMarkersDeseq2, sig_tissue_markergenes)
}

#We'll save the marker gene signatuers
save(mySignatures_RegionMarkersDeseq2, file = 'Signature_repository/RegionMarkersDeseq2.bin')

