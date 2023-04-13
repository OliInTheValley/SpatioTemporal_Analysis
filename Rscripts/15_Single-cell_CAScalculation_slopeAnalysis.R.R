myPaths <- .libPaths()
myPaths <- c(myPaths, '/home/olihahn/R/x86_64-pc-linux-gnu-library/3.6')
.libPaths(myPaths)

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(VISION)

###Having established that the CAS is changing with a region-specific magnitude in the bulkseq dataset, we now want to utilize single-cell/nuc datasets to find out
#which cell type(s) might be the one expressing the CAS and could therefore be 'responsible' for the increase on the bulk level
###To this wen, we will now load the signatures and calucalte scores for each cell in the single-cell transcriptome, using the signature library we generated before
#We will demonstrate this here using the hippocampus Nuc-seq dataset

seurat_object <- readRDS('R_objects/Nucseq_Hip_aging.rds')
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

scRNA_Aging_Scoretable <- seurat_object@meta.data
#We'll save the metadata of the object for later statistical analysis of the score cahgne
save(scRNA_Aging_Scoretable, file = 'R_objects/Nucseq_Hip_AgingDEGs_Scoretable.bin')

#Having mapped the scores for CAS back to the seurat object, we dan now visualize these using Seurat's standard functions
#This creates Figure 3E and Figure 3F, plotting the CAS directly in spatial location
VlnPlot(seurat_object, 'CommonAgingSignature', group.by = 'age', split.by = 'cell_type') + facet_grid(~split)


####We will make a copy of the original score table. Then we duplicate the values from the CommonAgingSignatue (CAS)
#To keep it simple, we'll just assign this column the label 'score'
score_table <- scRNA_Aging_Scoretable 
score_table$score <- score_table[,"CommonAgingSignature"]
#First we will have to make sure that 'age'  is set as a numberic vector. These get can messed up with all the different data transformations
score_table$age <-  as.numeric(plyr::mapvalues(score_table$age, from = c("Young", 'Old'), to = c(3,21)))

#Then we plot the scores for each age/cell type combination via violion plots
#This creates Figure 4B
(myplot <- ggplot(score_table, aes(x=age, y=score)) + geom_point(color='#bbbdbf', position = 'jitter') +geom_violin(notch = F, outlier.colour = NA,  mapping = aes(fill=as.factor(age))) +
    geom_smooth(data = score_table, mapping = aes(x=as.numeric(age), y=score), color='red', se = F, method = 'lm')  + 
    facet_grid(~cell_type) 
)

#Now we'll use the same stastitcal framework used for the bulkseq data to to properly qunatify and test the differences in CAS slopes
#First, we'll set the 'cell_type' informatino to a new, character-type column, as the factors can confuse the linear models
score_table$cell_type <- as.character(score_table$cell_type)
#First, we will calculate a linear model that allows for an interaction between age and cell_type (denoted with an)
m.interaction <- lm(score ~ age*cell_type, data = score_table)
#We are going to setup a dataframe that contains the intercepts for each cell_type for later inspection of the offset
coefficient_dataframe <- as.data.frame(m.interaction$coefficients)
colnames(coefficient_dataframe) <- 'coefficients'
coefficient_dataframe <- coefficient_dataframe[grep('^cell_type', row.names(coefficient_dataframe), value = T),,drop=F]
coefficient_dataframe$cell_type <- gsub('cell_type', '',row.names(coefficient_dataframe))
colnames(coefficient_dataframe)[1] <- 'offset'
#This should look like this:
#                                     offset              cell_type
# cell_typeBEC                     0.04650473                    BEC
# cell_typeChorPlexusEpithelial    0.04280939   ChorPlexusEpithelial
# cell_typeDG_Granule_neuron      -0.08285688      DG_Granule_neuron

#Now we'll use the function lstrends (fromt he lsmeans package) to quantify slopes with confidence intervals to allow cell type-to-cell type comparison
m.lst <- lstrends(m.interaction, "cell_type", var="age")
slope_dataframe <- as.data.frame(m.lst)
#Let's look at that output and store it in a separte dataframe to hold the information about the slopse
#Create pairwise comparisons between slopes and test using the paris function and store that result
#This lets us know which slopes are statistically significantly different from one another
slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))

#We will now plot the results from lstrends as bargraphs with confidence intervals
#We'll order the rank of the cell types in ascending order along their CAS slope estimate
slope_dataframe$cell_type <- factor(slope_dataframe$cell_type, levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = T),'cell_type']),
                                    ordered = T)
#We will color the bars accoring to the score incrase ('age.trend')
#This will create Figure 4C
(myplot <- ggplot(slope_dataframe, aes(x=cell_type, y=age.trend,fill=cell_type)) +  
    geom_bar(position=position_dodge(), stat="identity", mapping = aes(fill=cell_type)) + xlab("") + ylab("CAS slope")  +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    scale_y_continuous(expand = c(0,0)) + 
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1))
)


