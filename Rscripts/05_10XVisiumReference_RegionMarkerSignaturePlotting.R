
library(Seurat)
library(SeuratData)
setwd('/oak/stanford/scg/lab_twc/OHahn/RNA-seq/Documentation/2023_03_22_CAbulk_github/')

#We'll have to load the following datasets
#The coronal section needs to be downloaded from 10X directly: https://www.10xgenomics.com/resources/datasets/mouse-brain-section-coronal-1-standard-1-0-0
#Once loaded, the dataset should be stored in the input_data folder
brain_Coronal <- Load10X_Spatial("input_data/Visium_Reference/Mouse_coronal_unstained/filtered_feature_bc_matrix_Coronal_unstained/", filename = 'V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5')

#The saggital dataset is already available through the SeueratData package
InstallData('stxBrain.SeuratData')

#The saggital brain sections can also be optained dirclty throuhg the SeuratData  function
brain_Sag_Anterior <- LoadData("stxBrain.SeuratData", type = "anterior1")
brain_Sag_Posterior <- LoadData("stxBrain.SeuratData", type = "posterior1")

#First, we'll merge the two saggital sections into a single object
brain_Sag.merge <- merge(brain_Sag_Anterior, brain_Sag_Posterior)

#We'll filter spots that have fewer than 5 counts
brain_Sag.merge <- brain_Sag.merge[,brain_Sag.merge$nCount_Spatial >=5]
#To extract the sample label for each spot, we'll have to use bit of a combinatino of stringplit and lapply
brain_Sag.merge@meta.data$sample <-  unlist(lapply(strsplit(row.names(brain_Sag.merge@meta.data), split = "_"), function(x) x[[2]]))

#Relabel "samples"
brain_Sag.merge@meta.data$sampleID <-  plyr::mapvalues(x = brain_Sag.merge@meta.data$sample,
                                                       from = c('1',
                                                                '2'),
                                                       to = c('anterior',
                                                              'posterior')
)

# Now we'll integrate the transcriptomes to accomplish a better harmonization across the two independent capture areas
###Integration - Sagittal anterior/posterior
DefaultAssay(brain_Sag.merge) <- 'Spatial'
#We utilize Seurat's standard workflow for integration using SCT. First we'll split the object back into its original samples and perform SCTransform on these
object_splitlist <- SplitObject(brain_Sag.merge, split.by = "sample")
for (i in names(object_splitlist)) {
  object_splitlist[[i]] <- SCTransform(object_splitlist[[i]], verbose = T, assay = 'Spatial')
}
#We utilitze the selctIntegrationFeatures (with 2000 genes) and PrepSCTIntegration with default settings
Integration.features <- SelectIntegrationFeatures(object.list = object_splitlist, nfeatures = 2000)
object_splitlist <- PrepSCTIntegration(object.list = object_splitlist, anchor.features = Integration.features, verbose = T)
#the final Integration is then conducted
integration.anchors <- FindIntegrationAnchors(object.list = object_splitlist, normalization.method = "SCT",
                                                  anchor.features = Integration.features, verbose = T)
#Having identified the anchors, we can now run the integration
brain_Sag.merge.integrated <- IntegrateData(anchorset =integration.anchors, normalization.method = "SCT",
)
#The integration process creates some arbitrary, empty slots in the images tab - we'll remove those
brain_Sag.merge.integrated@images[-which(names(brain_Sag.merge.integrated@images) %in% c('anterior1', 'posterior1.1'))] <- NULL

#Now we will run the dimensionality reduction and clustering on the integrated data slot
DefaultAssay(brain_Sag.merge.integrated)
brain_Sag.merge.integrated <- RunPCA(object = brain_Sag.merge.integrated, verbose = T)
brain_Sag.merge.integrated <- FindNeighbors(brain_Sag.merge.integrated, dims = 1:30)
brain_Sag.merge.integrated <- FindClusters(brain_Sag.merge.integrated, resolution = 0.8)
brain_Sag.merge.integrated <- RunTSNE(brain_Sag.merge.integrated, dims = 1:30, verbose = T)
brain_Sag.merge.integrated <- RunUMAP(brain_Sag.merge.integrated, dims = 1:30, verbose = T)

#Let's check the outcome
#Let's look at the result - first just the original HE image
SpatialDimPlot(brain_Sag.merge.integrated, pt.size.factor = 0)
#And now we'll place transcriptome spots on top
SpatialDimPlot(brain_Sag.merge.integrated)
#We can also look just at transcriptome as a umap and visualize if and how similar anterior/posterior transcriptomes are
DimPlot(brain_Sag.merge.integrated, group.by = 'sample')

#Now we'll run SCTransformation as a way to get the expression data normalized in a way that we can analyze gene expression, calculate signature scores etc.

#Run SCT on the whole dataset to allow for differntial expression
#First we need to re-set the default assay ("integrated" is only sufficient for the dimensionality reduction)
DefaultAssay(brain_Sag.merge.integrated) <- "Spatial"
#Run SCT on the whole dataset
brain_Sag.merge.integrated <- SCTransform(brain_Sag.merge.integrated, assay = "Spatial", verbose = T)
#Set the default assay to SCT so it is the standard for all plotting/visualisations
DefaultAssay(brain_Sag.merge.integrated) <- 'SCT'
#We'll try out the Choroid Plexus marker gene Ttr as a test
SpatialFeaturePlot(brain_Sag.merge.integrated, 'Ttr')
#Let's save the integrated saggital dataset
saveRDS(brain_Sag.merge.integrated, file = 'R_objects/VisiumReference_Sagittal.rds')



###Now we run a similar analysis for the coronal section
DefaultAssay(brain_Coronal) <- "Spatial"

#filter for spots that don't have 0 counts
brain_Coronal <- brain_Coronal[,brain_Coronal$nCount_Spatial >=5]
#We only have one sample/slide here. So we can run SCTransform direclty and use this as input for dimensionality reduction
brain_Coronal <- SCTransform(brain_Coronal, assay = "Spatial", verbose = T)
#Make sure that SCTransfomred count matrix is set as default 
DefaultAssay(brain_Coronal) <- 'SCT'
#Now we'll run the dimensionality reduction and clustering
DefaultAssay(brain_Coronal)
brain_Coronal <- RunPCA(object = brain_Coronal, verbose = T)
brain_Coronal <- FindNeighbors(brain_Coronal, dims = 1:30)
brain_Coronal <- FindClusters(brain_Coronal, resolution = 0.8)
brain_Coronal <- RunUMAP(brain_Coronal, dims = 1:30, verbose = T)

#Let's look at the result - first just the original HE image
SpatialDimPlot(brain_Coronal, pt.size.factor = 0)
#And now we'll place transcriptome spots on top
SpatialDimPlot(brain_Coronal)
#We'll check out the same choroid plexus marker gene Ttr
SpatialFeaturePlot(brain_Coronal, 'Ttr')
saveRDS(brain_Coronal, file = 'R_objects/VisiumReference_Coronal.rds')

###Now we'll run the score analysis for the bulk marker genes that we discovered in 03_BulkSeq_DiagnosticUMAP_RegionMarkerAnalysis.R 
#for both the coronal and sagittal dataset
#First, we load the bulk region signature list
load('Signature_repository/RegionMarkersDeseq2.bin')

#First we will run the analysis for the coronal section
#We make sure SCT is set as default
DefaultAssay(brain_Coronal) <- 'SCT'
#We perpare the VISION object
Coronal_vision.obj <- Vision(brain_Coronal, signatures = mySignatures_RegionMarkersDeseq2, assay='SCT')
#We run the vision analysis
Coronal_vision.obj <- analyze(Coronal_vision.obj)

#The analysis has finished, and we'll remap now the scores for each of the signatures over to the metadata of the seurat object to keep them in one place
#We'll rotate over all the columns in the SigScore section of the vision object
for (signature_set in colnames(Coronal_vision.obj@SigScores)) {
  
  sigScores <- Coronal_vision.obj@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  #We use the spotID stored in the meta.data's rownames as an 'anchor' to map over the scores for the current signature
  brain_Coronal@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(brain_Coronal@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

## Now we'll repeat the same process for the sagittal, integrated dataset
#We make sure SCT is set as default
DefaultAssay(brain_Sag.merge.integrated) <- 'SCT'
#We perpare the VISION object
Sagittal_vision.obj <- Vision(brain_Sag.merge.integrated, signatures = mySignatures_RegionMarkersDeseq2, assay='SCT')
#We run the vision analysis
Sagittal_vision.obj <- analyze(Sagittal_vision.obj)

#The analysis has finished, and we'll remap now the scores for each of the signatures over to the metadata of the seurat object to keep them in one place
for (signature_set in colnames(Sagittal_vision.obj@SigScores)) {
  
  sigScores <- Sagittal_vision.obj@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  brain_Sag.merge.integrated@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(brain_Sag.merge.integrated@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

#Having mapped over the signatures, we can plot each signature similar to plotting a single gene, since seurat can access numeric values present in the meta data
#We'll look now how the signature of the white matter/corpus callosum ('c') look like in coronal data ....
SpatialFeaturePlot(brain_Coronal, 'CA1_cc_RegionMarkersDeseq2')
#... and in the sagittal dataset
SpatialFeaturePlot(brain_Sag.merge.integrated, 'CA1_cc_RegionMarkersDeseq2')
#This creates the plots for Figure S3A
#We'll save these two seurat objects
saveRDS(brain_Sag.merge.integrated, file = 'R_objects/VisiumReference_Sagittal.rds')
saveRDS(brain_Coronal, file = 'R_objects/VisiumReference_Coronal.rds')
