library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


###We have detected common and unique signatures of brain rejuvenation and aging We will now load the signatures and calculate scores for each rejuvenation intervention and region, 
#using the signature library we generated before and the integrated seurat object which we used for the quality control of the rejuvenation samples

seurat_object <- readRDS('R_objects/CA2_seurat.rds')
#We load the aging signature set, as well as both rejuvenation signature sets and construct a signature list out of these three
load('Signature_repository/AgingCohort_BulkDEG_signatures.bin')
load('Signature_repository/RejuvenationCohort_aDR_signatures.bin')
load('Signature_repository/RejuvenationCohort_YMP_signatures.bin')

signature_list <- c(mySignatures_CAS_AgingDEGs_perRegion, mySignatures_aDR_DEG, mySignatures_YMP_DEG)

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

BulkSeq_Rejuvenation_Scoretable <- seurat_object@meta.data
#We'll save the metadata of the object for later statistical analysis
save(BulkSeq_Rejuvenation_Scoretable, file = 'R_objects/BulkSeq_Rejuvenation_Scoretable.bin')





