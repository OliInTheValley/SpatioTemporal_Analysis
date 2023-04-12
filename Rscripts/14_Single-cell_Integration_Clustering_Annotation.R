
#Here, we are performing the single-cell/nuc RNA-seq analysis of the hippocampus. This consists of loading the pre-processed data that can be obtained from the GEO
#We load the individual samples and assemble them into a single object. Since these samples contain nuclei collections of males and femlaes,
#we'll use sex-specific transcripts to demultiplex these. Then we run a default integration, clustering and marker analysis, before annotating the cell clusters
#This is the defual workflow for the hippocampus and caudate putamen. Besides of the demultiplexing step, this workflow is also used for all other single-nuc/cell datasets utilized in this study.
#We will therefore only demonstrate this anaysis once using the hippocampus and recommend the user to adapt it to their dataset of interest

library(Seurat)
library(SeuratData)
setwd('/oak/stanford/scg/lab_twc/OHahn/RNA-seq/Documentation/2023_03_22_CAbulk_github/')

#We'll have to load the following datasets
#This works best if the user organises the per-sample package of barcode, matrix and feature files into sub-directories 
#that are labelled in some way that it represents region, age and replicate
#An example could look like this for GSM6537880 (Young hippocampus replicate 01)
#input_data/Nucseq_Hippocampus/HIP_Y_1/barcodes.tsv.gz
#input_data/Nucseq_Hippocampus/HIP_Y_1/features.tsv.gz
#input_data/Nucseq_Hippocampus/HIP_Y_1/matrix.mtx.gz

###
dataset_1_dir <- "input_data/Nucseq_Hippocampus/HIP_Y_1/"
data_young_01 <- Read10X(data.dir = dataset_1_dir)
data_young_01 <- CreateSeuratObject(counts = data_young_01, min.cells = 3, min.features = 40)
data_young_01$sampleID <- 'Hip_Young_R01'

dataset_2_dir <- "input_data/Nucseq_Hippocampus/HIP_Y_2/"
data_young_02 <- Read10X(data.dir = dataset_2_dir)
data_young_02 <- CreateSeuratObject(counts = data_young_02, min.cells = 3, min.features = 40)
data_young_02$sampleID <- 'Hip_Young_R01'

dataset_3_dir <- "input_data/Nucseq_Hippocampus/HIP_O_1/"
data_old_01 <- Read10X(data.dir = dataset_3_dir)
data_old_01 <- CreateSeuratObject(counts = data_old_01, min.cells = 3, min.features = 40)
data_old_01$sampleID <- 'Hip_Old_R01'

dataset_4_dir <- "input_data/Nucseq_Hippocampus/HIP_O_2/"
data_old_02 <- Read10X(data.dir = dataset_4_dir)
data_old_02 <- CreateSeuratObject(counts = data_old_02, min.cells = 3, min.features = 40)
data_old_02$sampleID <- 'Hip_Old_R02'

#We will put the individual samples together into a simple list
Nuc_Seq_list <- list()
Nuc_Seq_list[['Hip_Young_R01']] <- data_young_01
Nuc_Seq_list[['Hip_Young_R02']] <- data_young_02
Nuc_Seq_list[['Hip_Old_R01']] <- data_old_01
Nuc_Seq_list[['Hip_Old_R02']] <- data_old_02

#We will loop through the list and perform a few simple relabeling and filtering steps
#Notably, we remove the gene Malat1 from the count matrix. It's a highly expressed non-coding RNA that can cause issues during normalization
#If the user perfers to keep it in, just assign an empty vector to the 'genes_to_remove' object
genes_to_remove <- c("Malat1")

#Loop over the other datasets and place in a list
for (sample_to_analyze in names(Nuc_Seq_list) ) {
  print(sample_to_analyze)
  #We'll extract the respective single cell object
  sc_RNAseq_object <- Nuc_Seq_list[[sample_to_analyze]]
  #Here we remove any genes that we defined before starting the loop
  if (length(genes_to_remove) > 0) {
    #We extract the count matrix of the respective sample
    counts <- GetAssayData(sc_RNAseq_object, assay = "RNA")
    #subset the count matrix to all genes except those in the genes_to_remove vector
    counts <- counts[-(which(rownames(counts) %in% genes_to_remove)),]
    #Then subset the scRNA object to contain only genes that exist in the subsetted count matrix
    sc_RNAseq_object <- subset(sc_RNAseq_object,features = rownames(counts) )
  }
  #We re-create the sc object
  cleaned_matrix <- as.matrix(sc_RNAseq_object@assays$RNA@counts)
  sc_RNAseq_object <- CreateSeuratObject(counts = cleaned_matrix, min.cells = 3, min.features = 50)
  sc_RNAseq_object@meta.data$sample <- sample_to_analyze
  #we quantify the percent of reads mapping to mitochondrial and ribosomal genes. Mitochondrial reads are an indicator of damaged cells/nuclei
  sc_RNAseq_object[["percent.mt"]] <- PercentageFeatureSet(sc_RNAseq_object, pattern = "mt-")
  sc_RNAseq_object[["percent.rpl"]] <- PercentageFeatureSet(sc_RNAseq_object, pattern = "rp-")
  #We filter for cells/nuclei that express more than 400 genes but less than 7000, and retain nuclei with less than 5% mitochondrial reads
  sc_RNAseq_object <- subset(sc_RNAseq_object, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 5)
  #We'll assign the sample name into a new column 'sampleID'
  sc_RNAseq_object$sampleID <- sample_to_analyze
  
  ####Now, we'll demultiplex the single-cell object based on expression of sex-specific genes. 
  #This step is only necessary for the caudate putamen and hippocampus aging nuc-seq data
  #To this end, we perform a standard normalization 
  #sc_RNAseq_object <- NormalizeData(sc_RNAseq_object)
  #sc_RNAseq_object <- ScaleData(sc_RNAseq_object, features = rownames(sc_RNAseq_object))
  
  #Calculate the percentage of female-specific genes
  #sc_RNAseq_object[["percent.female"]] <- PercentageFeatureSet(sc_RNAseq_object, pattern = "Xist|Tsix")
  ##Calculate the percentage of male-specific genes
  sc_RNAseq_object[["percent.male"]] <- PercentageFeatureSet(sc_RNAseq_object, pattern = "Ddx3y|Eif2s3y|Uty|Kdm5d")
  #Calculate the log2Foldchange between the % of female to male reads
  #sc_RNAseq_object[["log2FtM"]] <- log2(sc_RNAseq_object[["percent.female"]]/sc_RNAseq_object[["percent.male"]])
  #Now we assign each cell if it's derived from the male or female mouse based on a +1/-1 log2FC cutoff. Everything in between is labelled 'unclear'
  #sc_RNAseq_object[['type']] <- ifelse(is.na(sc_RNAseq_object@meta.data$log2FtM), 'unclear', ifelse(sc_RNAseq_object@meta.data$log2FtM > 1, 'female', ifelse(sc_RNAseq_object@meta.data$log2FtM < -1, 'male', 'unclear')))
  #We create a little diagnostic plot
  #myplot <- ggplot(data=sc_RNAseq_object@meta.data, mapping=aes(x=nCount_RNA, y=log2FtM, color=type)) + geom_point()
  #print(myplot)
  #print(table(sc_RNAseq_object[['type']]))
  #We then filter out all 'unclear' cells
  #sc_RNAseq_object <- subset(sc_RNAseq_object, subset= type!='unclear')
  sc_RNAseq_object[['sampleID_sex']] <- paste(sc_RNAseq_object$sampleID, sc_RNAseq_object$type, sep = '_')
  
  #Now we assign the object back into the main list 
  Nuc_Seq_list[[sample_to_analyze]] <- sc_RNAseq_object
  rm(sc_RNAseq_object)
  
}

#Having prepared all samples individually, we can begin the the integration steps
#First, we'll merge all the individual samples into a single object
scRNA.merge <- Reduce(function(x,y) merge(x,y) , Nuc_Seq_list) 
rm(Nuc_Seq_list)

# Now we'll integrate the transcriptomes to accomplish a better harmonization across all individual samples
###Integration
DefaultAssay(scRNA.merge) <- 'RNA'
#We utilize Seurat's standard workflow for integration using SCT. First we'll split the object back into its original samples and perform SCTransform on these
object_splitlist <- SplitObject(scRNA.merge, split.by = "sampleID")
for (i in names(object_splitlist)) {
  print(i)
  object_splitlist[[i]] <- SCTransform(object_splitlist[[i]], verbose = T, assay = 'RNA', vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), )
}
#We utilitze the selctIntegrationFeatures (with 500 anchor genes) and PrepSCTIntegration with default settings
Integration.features <- SelectIntegrationFeatures(object.list = object_splitlist, nfeatures = 500)
object_splitlist <- PrepSCTIntegration(object.list = object_splitlist, anchor.features = Integration.features, verbose = T)
#the final Integration is then conducted
integration.anchors <- FindIntegrationAnchors(object.list = object_splitlist, normalization.method = "SCT",
                                              anchor.features = Integration.features, verbose = T)
#Having identified the anchors, we can now run the integration
scRNA.merge.integrated <- IntegrateData(anchorset =integration.anchors, normalization.method = "SCT",
)

#Before we move forward, we do a bit of relabeling and add informatino about the age groups
scRNA.merge.integrated@meta.data$age <-  plyr::mapvalues(x = scRNA.merge.integrated@meta.data$sampleID,
                                                         from = c('Hip_Young_R01',
                                                                  'Hip_Young_R02',
                                                                  'Hip_Old_R01',
                                                                  'Hip_Old_R02'
                                                                  ),
                                                         to = c('Young',
                                                                'Young',
                                                                'Old',
                                                                'Old')
)
#it is a bit easier to work with for now to use discrete values for the age groups instead of using age as a continious variable
scRNA.merge.integrated@meta.data$age <- factor(scRNA.merge.integrated@meta.data$age, 
                                               levels = c( 
                                                 'Young',
                                                 'Old'
                                               ), 
                                               ordered = T)


#we'll also extract the replicate information
scRNA.merge.integrated@meta.data$replicate <- unlist(lapply(strsplit((scRNA.merge.integrated@meta.data$sampleID), split = "_"), function(x) x[[3]]))
#and the cell BC
scRNA.merge.integrated@meta.data$index <- unlist(lapply(strsplit(row.names(scRNA.merge.integrated@meta.data), split = "-"), function(x) x[[1]]))
#Finlally, we'll create a composite cell index based on sampleID and cell BC
scRNA.merge.integrated@meta.data$cell_index <- paste(scRNA.merge.integrated@meta.data$tissue, scRNA.merge.integrated@meta.data$age, scRNA.merge.integrated@meta.data$replicate, scRNA.merge.integrated@meta.data$index, sep = '_')

#We'll also reassing some information that we already have to make the data better match to the metadata of bulk and spatial data
scRNA.merge.integrated@meta.data$tissue <- "Hip"
scRNA.merge.integrated@meta.data$sex <- scRNA.merge.integrated@meta.data$type


#Now we will run the dimensionality reduction and clustering on the integrated data slot
DefaultAssay(scRNA.merge.integrated)
scRNA.merge.integrated <- RunPCA(object = scRNA.merge.integrated, verbose = T)
scRNA.merge.integrated <- FindNeighbors(scRNA.merge.integrated, dims = 1:12)
scRNA.merge.integrated <- FindClusters(scRNA.merge.integrated, resolution = 0.4)
scRNA.merge.integrated <- RunUMAP(scRNA.merge.integrated, dims = 1:12, verbose = T)

#Let's check the outcome


#We can also look just at transcriptome as a umap and visualize if and how similar anterior/posterior transcriptomes are
DimPlot(scRNA.merge.integrated, label = T)
DimPlot(scRNA.merge.integrated, group.by = 'sampleID')

#Now we'll run SCTransformation as a way to get the expression data normalized in a way that we can analyze gene expression, calculate signature scores etc.

#First we need to re-set the default assay ("integrated" is only sufficient for the dimensionality reduction)
DefaultAssay(scRNA.merge.integrated) <- "RNA"
#Run SCT on the whole dataset
scRNA.merge.integrated <- SCTransform(scRNA.merge.integrated, assay = "RNA", verbose = T)
#In addtion, we'll also run the standard data normalization and scaling
DefaultAssay(scRNA.merge.integrated) <- "RNA"

scRNA.merge.integrated <- NormalizeData(scRNA.merge.integrated)
scRNA.merge.integrated <- FindVariableFeatures(scRNA.merge.integrated, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA.merge.integrated)
scRNA.merge.integrated <- ScaleData(scRNA.merge.integrated, features = all.genes)

#Let's save the integrated dataset
saveRDS(scRNA.merge.integrated, file = 'R_objects/Nucseq_Hip_aging.rds')


### Next, we need to annotate the clusters. For that we perform marker analysis for each cluster
Idents(scRNA.merge.integrated) <- 'seurat_clusters'
MarkerList_to_annotate <- FindAllMarkers(object = scRNA.merge.integrated, only.pos = TRUE, min.pct = 0.15, 
                                         thresh.use = 0.15, verbose = T, assay = 'SCT')


####We'll go now through the marker genes and identify which clusters represent which cell types based on these

####Having identified the the most likely cell types, we will now annotate the clusters with region labels
####We have identified 24 clusters - there are several cluster representing doublets, marked by expressing cell type markers of differing cell types
#We also found one small populations which could represent a special type of astrocytes or astrocyte-astrocyte doublets.

current.cluster.ids <- c(0, 1, 2, 3,
                         4, 5, 6, 7,
                         8, 9, 10, 11,
                         12, 13, 14, 15,
                         16, 17, 18, 19,
                         20, 21, 22, 23)

new.cluster.ids <- c('DG_Granule_neuron', 'DG_Granule_neuron', 'Astrocyte', 'DG_Granule_neuron',
                         'Mature_Oligodendrocyte', 'DG_Granule_neuron', 'Excitatory_neuron', 'Microglia',
                         'Interneuron', 'Excitatory_neuron', 'Doublet', 'Interneuron',
                         'Excitatory_neuron', 'Unkown', 'OPC', 'Interneuron',
                         'ChorPlexusEpithelial', 'BEC', 'Excitatory_neuron', 'Doublet',
                         'Doublet', 'Doublet', 'Doublet', 'Doublet')

#we use plyr's mapvalues function to re-label the clusters, 
scRNA.merge.integrated@meta.data$cell_type <- plyr::mapvalues(x = scRNA.merge.integrated@meta.data$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

#Let's set the default 'idents' to our cluster-level annotation
Idents(scRNA.merge.integrated) <- 'cell_type'
#Let's visualize the results as umap
#These function create the single-cell/nuc panesl for Figures 4,6, S9, S10, S11
DimPlot(scRNA.merge.integrated)
FeaturePlot(scRNA.merge.integrated, 'C4b', split.by = 'age')
VlnPlot(scRNA.merge.integrated, 'C4b', group.by = 'cell_type', split.by = 'age')

#We'll remove the doublets and unknown clusters
scRNA.merge.integrated <- subset(scRNA.merge.integrated, cell_type != 'Doublet')
scRNA.merge.integrated <- subset(scRNA.merge.integrated, cell_type != 'Unkown')
#We'll save the object
saveRDS(scRNA.merge.integrated, file = 'R_objects/Nucseq_Hip_aging.rds')


#Now we have everything to perform differential expression between young and old
#First we create a new categorial variable in the metadata, which combines the information of region and age
scRNA.merge.integrated@meta.data$cell_type_age <-
  paste(scRNA.merge.integrated@meta.data$cell_type,
        scRNA.merge.integrated@meta.data$age,
        sep = '-'
  )

#We set the default assay to 'spatial
DefaultAssay(scRNA.merge.integrated) <- 'RNA'
#and the default idents to our newly-set up cell_type_age 
Idents(scRNA.merge.integrated) <- 'cell_type_age'
#Test used is DESEe2
test_to_use <- 'MAST'
#We also define as the minimal percent a gene needs to be expressed as 10% in at least one of the compared spots (e.g. at least 10% in the transcriptome spots coming from a young cortex)
min_pct <- 0.05
logfc_threshold <- 0.2
assay_to_use <- 'RNA'
latent.vars = NULL
random.seed <- 42

#We'll setup a list that will store the results from all cell_types tested
results_list_scRNA.merge.integrated <- list()

#Now we'll loop over the two cell_types cortex and white matter, perform differential expression analysis and store the results
for (cell_type_to_test in unique(scRNA.merge.integrated$cell_type)) {
  #We'll set up a temporatory list to store the results tables
  results_list_cell_typeLevel <- list()
  print(as.character(cell_type_to_test))
  #We set the ages we want to compare
  cond1 <- "Old"
  cond2 <- "Young"
  #concatenate that with the currently analysed cell_type
  cell_type_cond1 <- paste(cell_type_to_test, cond1, sep = '-')
  cell_type_cond2 <- paste(cell_type_to_test, cond2, sep = '-')
  print(paste(cell_type_cond1, 'vs', cell_type_cond2))
  #Then we perform the differential expression analysis with the settings defined above
  diffmarkers_table <-
    FindMarkers(
      object = scRNA.merge.integrated,
      ident.1 = cell_type_cond1,
      ident.2 = cell_type_cond2,
      test.use = test_to_use,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      assay = assay_to_use,
      verbose = T,
      latent.vars = latent.vars
    )
  #When Deseq2 is set as test, the output is already corrected for multiple testing
  #We store the genes (currently only in rownames) to a seprate column
  diffmarkers_table$gene_symbol <- row.names(diffmarkers_table)
  #The default padjust method for MAST is Bonferroni correction, which is a valid albeit relatively harsh method for this many genes
  #We will therefore also perform the BH correction
  diffmarkers_table$p_adjBH <- p.adjust(diffmarkers_table$p_val, method = 'BH')
  #we write out the comparison so that we can access the results table easily later
  comparison <- paste(cond1, cond2, sep = '_vs_')
  #Then we store the table containing the differnetial expresison results in the resultslist
  results_list_cell_typeLevel[[comparison]]$resall <- diffmarkers_table
  #For easier inspection of the results, we also store a subset of the table, where we already selected for genes passing the significance threshold of 0.05
  resSig <- subset(diffmarkers_table, p_val_adj < 0.05)
  results_list_cell_typeLevel[[comparison]]$ressig <- resSig
  
  #We store this temporary resultslist in the main list we created before starting the loop and save this before starting with the next tissue/cell_type
  results_list_scRNA.merge.integrated[[as.character(cell_type_to_test)]] <- results_list_cell_typeLevel
  save(results_list_scRNA.merge.integrated, 
       file='R_objects/Nucseq_Hip_AgingDEGs_results.bin')
  
  
}

