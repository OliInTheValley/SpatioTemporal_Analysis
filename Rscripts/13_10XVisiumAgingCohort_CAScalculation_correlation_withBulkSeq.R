
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

###Having established that the CAS is changing with a region-specific magnitude in the bulkseq dataset, we want to analyze here if and to what extend this observation can be replicated with spatial transcitomics
###To this wen, we will now load the signatures and calculate scores for each spot in the spatial transcriptome, using the signature library we generated before

seurat_object <- readRDS('R_objects/Visium_AgingCohort_Coronal.rds')
load('Signature_repository/AgingCohort_BulkDEG_signatures.bin')
signature_list <- mySignatures_CAS_AgingDEGs_perRegion

#Since this is bulk data we want to make sure the default assay should be set to SCT
DefaultAssay(seurat_object) <- 'SCT'
#We perpare the VISION object
BulkSeq_vision.obj <- Vision(seurat_object, signatures = signature_list, assay='SCT')
#We run the vision analysis
BulkSeq_vision.obj <- analyze(BulkSeq_vision.obj)

#The analysis has finished, and we'll remap now the scores for each of the signatures over to the metadata of the seurat object to keep them in one place
for (signature_set in colnames(BulkSeq_vision.obj@SigScores)) {
  
  sigScores <- BulkSeq_vision.obj@SigScores[,signature_set]
  sigRanks <- rank(sigScores)
  seurat_object@meta.data[,signature_set] <- as.numeric(plyr::mapvalues(x = row.names(seurat_object@meta.data), from = names(sigScores), to = as.numeric(sigScores)))
  
}

Visium_Aging_Scoretable <- seurat_object@meta.data
#We'll save the metadata of the object for later statistical analysis of the score cahgne
save(Visium_Aging_Scoretable, file = 'R_objects/Visium_Aging_Scoretable.bin')

#Having mapped the scores for CAS back to the seurat object, we dan now visualize these using Seurat's standard functions
#This creates Figure 3E and Figure 3F, plotting the CAS directly in spatial location
SpatialFeaturePlot(seurat_object, 'CommonAgingSignature', min.cutoff = 0, max.cutoff = 1)
VlnPlot(seurat_object, 'CommonAgingSignature', group.by = 'age', split.by = 'regionLevel') + facet_grid(~split)


####We will make a copy of the original score table. Then we duplicate the values from the CommonAgingSignatue (CAS)
#To keep it simple, we'll just assign this column the label 'score'
score_table <- Visium_Aging_Scoretable 
score_table$score <- score_table[,"CommonAgingSignature"]
#First we will have to make sure that 'age'  is set as a numberic vector. These get can messed up with all the different data transformations
score_table$age <-  as.numeric(gsub('M', "",as.character(score_table$age)))

#Then we plot the scores for each region along the age-axis. We also quanitfy the per-age group in a boxplot
#This creates Figure 3F
(myplot <- ggplot(score_table, aes(x=age, y=score)) + geom_point(color='#bbbdbf', position = 'jitter') +geom_violin(notch = F, outlier.colour = NA,  mapping = aes(fill=as.factor(age))) +
    geom_smooth(data = score_table, mapping = aes(x=as.numeric(age), y=score), color='red', se = F, method = 'lm')  + 
    facet_grid(~regionLevel) 
)

#Now we'll use the same stastitcal framework used for the bulkseq data to to properly qunatify and test the differences in CAS slopes
#First, we'll set the 'regionLevel' informatino to a new, character-type column, as the factors can confuse the linear models
score_table$tissue <- as.character(score_table$regionLevel)
#First, we will calculate a linear model that allows for an interaction between age and tissue (denoted with an)
m.interaction <- lm(score ~ age*tissue, data = score_table)
#We are going to setup a dataframe that contains the intercepts for each tissue for later inspection of the offset
coefficient_dataframe <- as.data.frame(m.interaction$coefficients)
colnames(coefficient_dataframe) <- 'coefficients'
coefficient_dataframe <- coefficient_dataframe[grep('^tissue', row.names(coefficient_dataframe), value = T),,drop=F]
coefficient_dataframe$tissue <- gsub('tissue', '',row.names(coefficient_dataframe))
colnames(coefficient_dataframe)[1] <- 'offset'
#This should look like this:
#                                          offset                        tissue
# tissueBasolateral amygdalar nucleus -0.0429308716 Basolateral amygdalar nucleus
# tissueChoroid Plexus                 0.2959964855                Choroid Plexus
# tissueCortex                        -0.0063221670                        Cortex
# tissueCortical subplate              0.0006455123             Cortical subplate
# tissueGlobus pallidus                0.1500761701               Globus pallidus


#Now we'll use the function lstrends (fromt he lsmeans package) to quantify slopes with confidence intervals to allow tissue-to-tissue comparison
m.lst <- lstrends(m.interaction, "tissue", var="age")
#Let's look at that output and store it in a separte dataframe to hold the information about the slopse
Visium_slope_dataframe <- as.data.frame(m.lst)

#We take the coeefiecnet dataframe that we gnereate above together with the slope dataframe
Visium_coefficient_dataframe <- merge.data.frame(slope_dataframe, coefficient_dataframe, by = 'tissue')

#Now we are interested in comparing the slopse as quantified here with Visium are matching those found in bulkseq
#First, we'll need to relable the 'tissue'/regions in the slope table to match the acronyms of the bulkdata
Visium_slope_dataframe$tissue_acronym <- plyr::mapvalues(Visium_slope_dataframe$tissue, 
                                                               from = c('Cortex', 'Hippocampus', 'Hypothalamus', 'Striatum', 'Thalamus', 'White matter', 'Choroid Plexus'),
                                                               to = c('cor', 'hi', 'hy', 'cp', 'th', 'cc', 'plx'))

#We only did this for the seven regions that we can quantify well in the coronal section of the Visium data, so we'll subset the coefficient dataframe to these
Visium_slope_dataframe <- Visium_slope_dataframe %>% filter(tissue_acronym %in% c('cor', 'hi', 'hy', 'cp', 'th', 'cc', 'plx'))
#So we load the slope table that we got out from bulk seq
load('R_objects/BulkSeq_CAS_slopes.bin')
Bulk_slope_dataframe <- slope_dataframe %>% filter(tissue %in% c('cor', 'hi', 'hy', 'cp', 'th', 'cc', 'plx'))
colnames(Bulk_slope_dataframe)[1] <- 'tissue_acronym'
#Now we'll merge the two dataframes
plot_df <- merge.data.frame(Bulk_slope_dataframe, Visium_slope_dataframe, by = 'tissue_acronym')

#We plot the output which creates Figure 3G
ggplot(plot_df, aes(x=age.trend.x, y=age.trend.y)) + geom_point() + geom_smooth(method = 'lm', se = F)

#We can also test the relationship
cor.test(plot_df$age.trend.x, plot_df$age.trend.y, method = 'spearman')
cor.test(plot_df$age.trend.x, plot_df$age.trend.y, method = 'pearson')
summary(lm(age.trend.y ~ age.trend.x, data = plot_df))
#We identify a significant association across all three tests


