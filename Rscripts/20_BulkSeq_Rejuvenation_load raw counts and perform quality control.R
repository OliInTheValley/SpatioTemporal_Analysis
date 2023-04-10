rm(list=ls())


library(readr)
library(DESeq2)

###Start
#Load dataset from raw counttable
library(readr)
BulkSeq_Rejuvenation_Counttable <- read_delim("input_data/BulkSeq_Rejuvenation/BulkSeq_Rejuvenation_Counttable.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)
BulkSeq_Rejuvenation_Counttable <- as.data.frame(BulkSeq_Rejuvenation_Counttable)
row.names(BulkSeq_Rejuvenation_Counttable) <- BulkSeq_Rejuvenation_Counttable[,1] #turn gene names into rownames of the table
BulkSeq_Rejuvenation_Counttable[1] <- NULL

#Load meta data and create an experimentDesign frame that fits Deseq2's design
experimentDesign  <- data.frame(row.names = colnames(BulkSeq_Rejuvenation_Counttable)) #Create dataframe with sampleIDs as rownames

BulkSeq_Aging_metadata <- read_delim("input_data/BulkSeq_Rejuvenation/BulkSeq_Rejuvenation_metadata.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
BulkSeq_Aging_metadata <- as.data.frame(BulkSeq_Aging_metadata)

#Transfer information from the meta data table into a Deseq2-compatible experimentDesign structure
for (meta_column in colnames(BulkSeq_Aging_metadata)[-which(colnames(BulkSeq_Aging_metadata)=='sampleID')]) {
  experimentDesign[,meta_column] <- plyr::mapvalues(x = row.names(experimentDesign), from = BulkSeq_Aging_metadata[,'sampleID'], to = BulkSeq_Aging_metadata[,meta_column])
}

(experimentDesign <- as.data.frame(unclass(experimentDesign), row.names = row.names(experimentDesign)))

#For simplicity, we will create a column called 'tissue' which contains acronyms for the verbose tissue names currently provided in the metadata

long_labels <- c("Corpus callosum",
                 "Cerebellum",
                 "Motor cortex",
                 "Caudate putamen",
                 "Entorhinal cortex",
                 "Hippocampus (posterior)",
                 "Hippocampus (anterior)", 
                 "Hypothalamus", 
                 "Medulla",
                 "Olfactory bulb", 
                 "Choroid Plexus",
                 "Pons",
                 "Subventricular zone",
                 "Thalamus",
                 "Visual Cortex",
                 "Liver")

acronym_labels <- c("cc",
                    "cer",
                    "cor",
                    "cp",
                    "ent",
                    "hi2",
                    "hi", 
                    "hy", 
                    "med",
                    "olf", 
                    "plx",
                    "pon",
                    "svz",
                    "th",
                    "vis",
                    "LL")
#Relabel 
experimentDesign$tissue <- plyr::mapvalues(x = experimentDesign$tissueLong, from = long_labels, to = acronym_labels)


#Create Deseq2 object
#This will be the main processed data object to work from. We use the acronym 'CA' for "cerebrum aevum" (latin for aging brain)
#to label the samples from the rejuvenation bulkseq dataset, we use numer '2', hence CA2
dds_CA2 <- DESeqDataSetFromMatrix(countData=BulkSeq_Rejuvenation_Counttable, colData=experimentDesign, design = ~treatment + tissue ) 
#we first run the analysis using a simple model without an interaction term
design(dds_CA2) <-  ~treatment + tissue

#First we'll inspect if the bulkseq data generated here captured the same regions as in the aging bulkseq dataset

#We create a Seurat object using the counts from the 'All' Desea2 object list (that's where we stored data from all the tissues)
CA2_seurat <- CreateSeuratObject(counts(dds_CA2))
#We'll also carry over some of the metadata
CA2_seurat@meta.data[,4:10] <- as.data.frame(colData(dds_CA2))[c('tissueLong', 'sex', 'age', 'mouseID', 'tissue', 'experiment_group', 'treatment')]
#The dataset currently contians data from liver tissue of aDR/AL mice. We won't need that for now so we'll filter it out
CA2_seurat <- subset(CA2_seurat, subset= tissue != "LL")

#We perform standard Log-normalization
CA2_seurat <- NormalizeData(CA2_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
#We define the variable genes (2000 should do)
CA2_seurat <- FindVariableFeatures(CA2_seurat, selection.method = "vst", nfeatures = 2000)

#Now we'll scale the data using all genes in the count matrix as input
CA2_seurat <- ScaleData(CA2_seurat, features = rownames(CA2_seurat), verbose = T)

#Now we'll run the PCA and visualize the explained variance per PC
CA2_seurat <- RunPCA(CA2_seurat)
ElbowPlot(CA2_seurat, ndims = 40)

#Now we'll run the graph-based clustering across the first 40 PCs
CA2_seurat <- FindNeighbors(CA2_seurat, dims = 1:40)
CA2_seurat <- FindClusters(CA2_seurat, resolution = 0.8)
#For visualization, we'll also run a UMAP across the same parameters
CA2_seurat <- RunUMAP(CA2_seurat, dims = 1:40, verbose = T)
#We'll set the default idents to be the same as in the 'tissue' column
Idents(CA2_seurat) <- 'tissue'
#We will need this later for score calcualtion, so we save the object
saveRDS(CA2_seurat, file='R_objects/CA2_seurat.rds')

#Now we load the seurat object from the aging cohort to enable inspection about the consistency of tissue isolation
CA1_seurat <- readRDS('R_objects/CA1_seurat.rds')
CA1_seurat$experiment_group <- 'WT Aging'

#And now we'll merge both datasets
CAcombined <- merge(CA1_seurat, CA2_seurat)

#We perform standard Log-normalization
CAcombined <- NormalizeData(CAcombined, normalization.method = "LogNormalize", scale.factor = 10000)
#We define the variable genes (2000 should do)
CAcombined <- FindVariableFeatures(CAcombined, selection.method = "vst", nfeatures = 2000)

#We'll inspect the most variable features manually
#first we'll define the top 10 variable features
top10 <- head(VariableFeatures(CAcombined), 10)
#And then use Seurat variablefeatureplot function
plot1 <- VariableFeaturePlot(CAcombined)
#Let's label these points using the labelpoints function
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Now we'll scale the data using all genes in the count matrix as input
CAcombined <- ScaleData(CAcombined, features = rownames(CAcombined), verbose = T)

#Now we'll run the PCA and visualize the explained variance per PC
CAcombined <- RunPCA(CAcombined)
ElbowPlot(CAcombined, ndims = 40)

#Now we'll run the graph-based clustering across the first 40 PCs
CAcombined <- FindNeighbors(CAcombined, dims = 1:40)
CAcombined <- FindClusters(CAcombined, resolution = 0.8)
#For visualization, we'll also run a UMAP across the same parameters
CAcombined <- RunUMAP(CAcombined, dims = 1:40, verbose = T)
#We'll set the default idents to be the same as in the 'tissue' column
Idents(CAcombined) <- 'tissue'

#Let's visualise the outcome 
#####And plot Figure S15GC
DimPlot(CAcombined, group.by = 'tissue', label = T)

###This looks good. So we will now split the Deseq2 object into two dataset, one for the young mouse plasma experiments and one for the dietary interventions
#First we create the Deseq2 object contatining data for the DR animals
#We will have to remove some empty factor levels to avoid problems with Deseq2 down the road
dds_CA2_DR <- dds_CA2[, dds_CA2$experiment_group == "Dietary Intervention"]
dds_CA2_DR$experiment_group <- droplevels(dds_CA2_DR$experiment_group)
dds_CA2_DR$treatment <- droplevels(dds_CA2_DR$treatment)
dds_CA2_DR$tissue <- droplevels(dds_CA2_DR$tissue)
#Run Deseq2
#We'll save this object for now
save(dds_CA2_DR, file='R_objects/dds_BulkSeq_DR.bin')


#Now we'll do the same for the plasma-treated mice
#We will have to remove some empty factor levels to avoid problems with Deseq2 down the road
dds_CA2_YMP <- dds_CA2[, dds_CA2$experiment_group == "Young Plasma Injection"]
dds_CA2_YMP$experiment_group <- droplevels(dds_CA2_YMP$experiment_group)
dds_CA2_YMP$treatment <- droplevels(dds_CA2_YMP$treatment)
dds_CA2_YMP$tissue <- droplevels(dds_CA2_YMP$tissue)
#Run Deseq2
save(dds_CA2_YMP, file='R_objects/dds_BulkSeq_YMP.bin')

