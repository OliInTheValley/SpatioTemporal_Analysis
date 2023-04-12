#Here, we are runing the base analysis of the Visium timecourse of Figure 3
#We'll assmeble the processed data that can be optained from GEO, integrate the data, before running clustering, marker anlaysis and region annotation
library(Seurat)
library(SeuratData)
setwd('/oak/stanford/scg/lab_twc/OHahn/RNA-seq/Documentation/2023_03_22_CAbulk_github/')

#We'll have to load the following datasets
###
dataset_1_dir <- "input_data/Visium_AgingCohort/Young_R01_S1/outs/"
dataset_1_filename <- "filtered_feature_bc_matrix.h5" 
data_young_01 <- Load10X_Spatial(data.dir = dataset_1_dir, filename = dataset_1_filename)
data_young_01$sampleID <- 'Visium_Young_R01_S1'

dataset_2_dir <- "input_data/Visium_AgingCohort/Mid_R01_S1/outs/"
dataset_2_filename <- "filtered_feature_bc_matrix.h5" 
data_mid_01 <- Load10X_Spatial(data.dir = dataset_2_dir, filename = dataset_2_filename)
data_mid_01$sampleID <- 'Visium_Mid_R01_S1'

dataset_3_dir <- "input_data/Visium_AgingCohort/Old_R01_S1//outs/"
dataset_3_filename <- "filtered_feature_bc_matrix.h5" 
data_old_01 <- Load10X_Spatial(data.dir = dataset_3_dir, filename = dataset_3_filename)
data_old_01$sampleID <- 'Visium_Old_R01_S1'

dataset_4_dir <- "input_data/Visium_AgingCohort/Young_R02_S1/outs/"
dataset_4_filename <- "filtered_feature_bc_matrix.h5" 
data_young_02 <- Load10X_Spatial(data.dir = dataset_4_dir, filename = dataset_4_filename)
data_young_02$sampleID <- 'Visium_Young_R02_S1'

dataset_5_dir <- "input_data/Visium_AgingCohort/Mid_R02_S2/outs/"
dataset_5_filename <- "filtered_feature_bc_matrix.h5" 
data_mid_02 <- Load10X_Spatial(data.dir = dataset_5_dir, filename = dataset_5_filename)
data_mid_02$sampleID <- 'Visium_Mid_R02_S2'

dataset_6_dir <- "input_data/Visium_AgingCohort/Old_R02_S2/outs/"
dataset_6_filename <- "filtered_feature_bc_matrix.h5" 
data_old_02 <- Load10X_Spatial(data.dir = dataset_6_dir, filename = dataset_6_filename)
data_old_02$sampleID <- 'Visium_Old_R02_S2'

#We will put the individual samples together into a simple list
Spatial_Seq_list <- list()
Spatial_Seq_list[['Visium_Young_R01_S1']] <- data_young_01
Spatial_Seq_list[['Visium_Mid_R01_S1']] <- data_mid_01
Spatial_Seq_list[['Visium_Old_R01_S1']] <- data_old_01
Spatial_Seq_list[['Visium_Young_R02_S1']] <- data_young_02
Spatial_Seq_list[['Visium_Mid_R02_S2']] <- data_mid_02
Spatial_Seq_list[['Visium_Old_R02_S2']] <- data_old_02

#First, we'll merge all the individual samples into a single object
Visium_brain_aging.merge <- Reduce(function(x,y) merge(x,y) , Spatial_Seq_list) 

#We'll filter spots that have fewer than 5 counts
Visium_brain_aging.merge <- Visium_brain_aging.merge[,Visium_brain_aging.merge$nCount_Spatial >=5]

# Now we'll integrate the transcriptomes to accomplish a better harmonization across the two independent capture areas
###Integration 
DefaultAssay(Visium_brain_aging.merge) <- 'Spatial'
#We utilize Seurat's standard workflow for integration using SCT. First we'll split the object back into its original samples and perform SCTransform on these
object_splitlist <- SplitObject(Visium_brain_aging.merge, split.by = "sampleID")
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
Visium_brain_aging.merge.integrated <- IntegrateData(anchorset =integration.anchors, normalization.method = "SCT",
)
#The integration process creates some arbitrary, empty slots in the images tab - we'll remove those
Visium_brain_aging.merge.integrated@images[-which(names(Visium_brain_aging.merge.integrated@images) %in% c('slice1', 'slice1.1.1', 'slice1.2.2', 'slice1.3.3', 'slice1.4.4', 'slice1.5.5', 'slice1.6.6'))] <- NULL

#Before we move forward, we do a bit of relabeling and add informatino about the age groups
Visium_brain_aging.merge.integrated@meta.data$age <-  plyr::mapvalues(x = Visium_brain_aging.merge.integrated@meta.data$sampleID,
                                                         from = c('Visium_Young_R01_S1',
                                                                  'Visium_Mid_R01_S1',
                                                                  'Visium_Old_R01_S1',
                                                                  'Visium_Young_R02_S1',
                                                                  'Visium_Mid_R02_S2',
                                                                  'Visium_Old_R02_S2'
                                                         ),
                                                         to = c('M6',
                                                                'M18',
                                                                'M21',
                                                                'M6',
                                                                'M18',
                                                                'M21')
)
#it is a bit easier to work with for now to use discrete values for the age groups instead of using age as a continious variable
Visium_brain_aging.merge.integrated@meta.data$age <- factor(Visium_brain_aging.merge.integrated@meta.data$age, 
                                               levels = c( 
                                                 'M6',
                                                 'M18',
                                                 'M21'
                                               ), 
                                               ordered = T)


#Finally, we set an order for the samples so we always have samples of the same age group paired next to one another
Visium_brain_aging.merge.integrated@meta.data$sampleID <- factor(Visium_brain_aging.merge.integrated@meta.data$sampleID, 
                                                                 levels = c('Visium_Young_R01_S1',
                                                                            'Visium_Young_R02_S1',
                                                                            'Visium_Mid_R01_S1',
                                                                            'Visium_Mid_R02_S2',
                                                                            'Visium_Old_R01_S1',
                                                                            'Visium_Old_R02_S2'
                                                                 ), 
                                                                 ordered = T)

#Now we will run the dimensionality reduction and clustering on the integrated data slot
DefaultAssay(Visium_brain_aging.merge.integrated)
Visium_brain_aging.merge.integrated <- RunPCA(object = Visium_brain_aging.merge.integrated, verbose = T)
Visium_brain_aging.merge.integrated <- FindNeighbors(Visium_brain_aging.merge.integrated, dims = 1:30)
Visium_brain_aging.merge.integrated <- FindClusters(Visium_brain_aging.merge.integrated, resolution = 0.8)
Visium_brain_aging.merge.integrated <- RunUMAP(Visium_brain_aging.merge.integrated, dims = 1:30, verbose = T)

#Let's check the outcome
#Let's look at the result - first just the original HE image
SpatialDimPlot(Visium_brain_aging.merge.integrated, pt.size.factor = 0)
#And now we'll place transcriptome spots on top
SpatialDimPlot(Visium_brain_aging.merge.integrated)
#We can also look just at transcriptome as a umap and visualize if and how well the data of the individual samples overlapxf
DimPlot(Visium_brain_aging.merge.integrated)
DimPlot(Visium_brain_aging.merge.integrated, split.by = 'sampleID')

#Now we'll run SCTransformation as a way to get the expression data normalized in a way that we can analyze gene expression, calculate signature scores etc.

#Run SCT on the whole dataset to allow for differntial expression
#First we need to re-set the default assay ("integrated" is only sufficient for the dimensionality reduction)
DefaultAssay(Visium_brain_aging.merge.integrated) <- "Spatial"
#Run SCT on the whole dataset
Visium_brain_aging.merge.integrated <- SCTransform(Visium_brain_aging.merge.integrated, assay = "Spatial", verbose = T)
#Set the default assay to SCT so it is the standard for all plotting/visualisations
DefaultAssay(Visium_brain_aging.merge.integrated) <- 'SCT'
#We'll try out the Choroid Plexus marker gene Ttr as a test
SpatialFeaturePlot(Visium_brain_aging.merge.integrated, 'Ttr')
#We'll try out if one of the top common aging score genes, C4b, exhibits a differnt expression pattern with age
SpatialFeaturePlot(Visium_brain_aging.merge.integrated, 'C4b')
#Let's save the integrated  dataset
saveRDS(Visium_brain_aging.merge.integrated, file = 'R_objects/Visium_AgingCohort_Coronal.rds')


### Next, we need to annotate the clusters. For that we perform marker analysis for each cluster
Idents(Visium_brain_aging.merge.integrated) <- 'seurat_clusters'
MarkerList_to_annotate <- FindAllMarkers(object = Visium_brain_aging.merge.integrated, only.pos = TRUE, min.pct = 0.1, 
                                         thresh.use = 0.1, verbose = T, assay = 'SCT')


####We'll go now through the marker genes and identify clusters based on these

####Having identified the the most likely cluster-region associations, we will now annotate the clusters with region labels
####We have identified 30 clusters
#We create a vector with our 'current' labels - simply the running nubmers of the seurat clusters
current.cluster.ids <- c(0, 1, 2, 3,
                         4, 5, 6, 7,
                         8, 9, 10, 11,
                         12, 13, 14, 15,
                         16, 17, 18, 19,
                         20, 21, 22, 23,
                         24, 25, 26, 27,
                         28, 29, 30)


#And we have set the labels of these clusters. each cluster is assigned one label
#Similar to our bulk Seq dataset, we have used acronyms without any whitespaces to facilitate downstream analyses

current.cluster.ids <- c(0, 1, 2, 3,
                         4, 5, 6, 7,
                         8, 9, 10, 11,
                         12, 13, 14, 15,
                         16, 17, 18, 19,
                         20, 21, 22, 23,
                         24, 25, 26, 27,
                         28, 29, 30)

new.cluster.ids_clusterLevel <- c('White matter', 'Thalamus 1', 'Hypothalamus', 'Striatum',
                          'Layer V', 'Thalamus 2', 'Globus pallidus', 'Layer II',
                          'Amygdala 1', 'Layer III', 'Cortical subplate', 'Layer IV',
                          'Tissueborder', 'Molecular layer', 'Stratum radiatum', 'Layer I',
                          'Layer VI', 'Retrospinal area 1', 'Thalamic reticular nucleus', 'Entorhinal area',
                          'Choroid Plexus', 'Basolateral amygdalar nucleus', 'Endopiriform nucleus', 'Thalamus 3',
                          'Ventricle', 'Cortical subplate', 'Retrospinal area 2', 'Dentate gyrus',
                          'Amygdala 2', 'CA3', 'CA1')

#In addition we have grouped the clusters into meaningful regions to allow comparison with the bulkSeq data

new.cluster.ids_regionLevel <- c('White matter', 'Thalamus', 'Hypothalamus', 'Striatum',
                     'Cortex', 'Thalamus', 'Globus pallidus', 'Cortex',
                     'Amygdala', 'Cortex', 'Cortical subplate', 'Cortex',
                     'Border', 'Hippocampus', 'Hippocampus', 'Cortex',
                     'Cortex', 'Cortex', 'Thalamic reticular nucleus', 'Cortex',
                     'Choroid Plexus', 'Basolateral amygdalar nucleus', 'Cortex', 'Thalamus',
                     'Ventricle', 'Cortical subplate', 'Cortex', 'Hippocampus',
                     'Amygdala', 'Hippocampus', 'Hippocampus')

#we use plyr's mapvalues function to re-label the clusters, 
Visium_brain_aging.merge.integrated@meta.data$clusterLevel <- plyr::mapvalues(x = Visium_brain_aging.merge.integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids_clusterLevel)
Visium_brain_aging.merge.integrated@meta.data$regionLevel <- plyr::mapvalues(x = Visium_brain_aging.merge.integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids_regionLevel)

#For easier visualization, we will also set both region annotations as vectors with ordered vectors - that way we have clusters that belong to the same region, organized next to each other

Visium_brain_aging.merge.integrated@meta.data$clusterLevel <- factor(Visium_brain_aging.merge.integrated@meta.data$clusterLevel, levels = c('Layer I', 'Layer II', 'Layer III', 'Layer IV', 'Layer V','Layer VI', 'Retrospinal area 1', 'Retrospinal area 2', 'Entorhinal area', 'Endopiriform nucleus',
                                                                                                                                            'Cortical subplate', 
                                                                                                                                            'Thalamus 1', 'Thalamus 2', 'Thalamus 3',
                                                                                                                                            'Thalamic reticular nucleus',
                                                                                                                                            'Striatum',
                                                                                                                                            'CA1', 'CA3', 'Dentate gyrus',  'Molecular layer', 'Stratum radiatum',
                                                                                                                                            'Hypothalamus',
                                                                                                                                            'Amygdala 1', 'Amygdala 2',
                                                                                                                                            'Basolateral amygdalar nucleus',
                                                                                                                                            'Globus pallidus',
                                                                                                                                            'White matter',
                                                                                                                                            'Choroid Plexus', 
                                                                                                                                            'Ventricle', 
                                                                                                                                            'Tissueborder'), 
                                                                     ordered = T)

Visium_brain_aging.merge.integrated@meta.data$regionLevel <- factor(Visium_brain_aging.merge.integrated@meta.data$regionLevel, levels = c('Cortex', 'Cortical subplate', 'Thalamus', 'Hippocampus',
                                                                                'Hypothalamus', 'Thalamic reticular nucleus','Striatum', 'Amygdala', 'Basolateral amygdalar nucleus', 'Globus pallidus', 'White matter',
                                                                                'Choroid Plexus', 'Ventricle', 'Border'), 
                                       ordered = T)



#Let's set the default 'idents' to our cluster-level annotation
Idents(Visium_brain_aging.merge.integrated) <- 'clusterLevel'
#Let's visualize the results as spatial map and as umap
SpatialDimPlot(Visium_brain_aging.merge.integrated)
DimPlot(Visium_brain_aging.merge.integrated)

#We'll do the same for the regionLevel annotation 
SpatialDimPlot(Visium_brain_aging.merge.integrated, group.by = 'regionLevel')
DimPlot(Visium_brain_aging.merge.integrated, group.by = 'regionLevel')

#It ooks like there is a cluster lining primarily the tissueborder (i.e. where the spots are on the very outside of the tissue)
#These are often contaminated with spurious RNA that floats around during the permeabilitzeino
SpatialDimPlot(Visium_brain_aging.merge.integrated, cells.highlight = CellsByIdentities(object = Visium_brain_aging.merge.integrated, idents = c('tissueborder')))
#We'll filter these spots out 
Visium_brain_aging.merge.integrated <- subset(Visium_brain_aging.merge.integrated, clusterLevel != 'Tissueborder')
#We'll save the object
saveRDS(Visium_brain_aging.merge.integrated, file = 'R_objects/Visium_AgingCohort_Coronal.rds')

#We create Figure S5 using the following commands
SpatialDimPlot(Visium_brain_aging.merge.integrated, group.by = 'regionLevel') #Visualize clustering
DimPlot(Visium_brain_aging.merge.integrated, group.by = 'regionLevel')
SpatialDimPlot(Visium_brain_aging.merge.integrated, group.by = 'regionLevel', split) #Visualize clustering
DimPlot(Visium_brain_aging.merge.integrated, group.by = 'regionLevel')



