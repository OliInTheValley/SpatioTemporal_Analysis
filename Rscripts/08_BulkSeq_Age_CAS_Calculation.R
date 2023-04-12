library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


###We have detected common and unique signatures of brain aging. We will now load the signatures and calculate scores for each, using the signature library we generated before
#and the integrated seurat object which we used for the quality control and region-marker analysis
seurat_object <- readRDS('R_objects/CA1_seurat.rds')
load('Signature_repository/AgingCohort_BulkDEG_signatures.bin')
signature_list <- mySignatures_CAS_AgingDEGs_perRegion

#Since this is bulk data we want to make sure the default assay should be set to RNA
DefaultAssay(seurat_object) <- 'RNA'
#We perpare the VISION object
BulkSeq_vision.obj <- Vision(seurat_object, signatures = signature_list, assay='RNA')
#We run the vision analysis
BulkSeq_vision.obj <- analyze(BulkSeq_vision.obj)

#The analysis has finished, and we'll remap now the scores for each of the signatures over to the metadata of the seurat object to keep them in one place
for (signature_set in colnames(BulkSeq_vision.obj@SigScores)) {
  
  sigScores <- BulkSeq_vision.obj@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  seurat_object@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(seurat_object@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

BulkSeq_Aging_Scoretable <- seurat_object@meta.data
#We'll save the metadata of the object for later statistical analysis of the score cahgne
save(BulkSeq_Aging_Scoretable, file = 'R_objects/BulkSeq_Aging_Scoretable.bin')





