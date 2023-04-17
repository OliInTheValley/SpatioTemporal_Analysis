##In the previous analyses, we have identified and applied gene signatures found across brain regions, that exhibit region-specific magnitudes in terms of how they change with age
##Before moving forward, we want to assess the possibliity of differences in cell composition would be a possible explanation for these differences
#That is, it could be that a different abundance of e.g. microglia at young could be the cause of why we observed stronger up-regulateion of the CAS
#To assess this, we took advantage of a pervioulsy-published STARmap dataset from Shi et al, 2022 (https://www.biorxiv.org/content/10.1101/2022.06.20.496914v1) by the lab of Xiao Wang.
#There, the authors quanitfied cell composition across a whole mouse brain using STARmap PLUS. We obtained the meta-data of their processed datasets, which we downloaded here:
#https://singlecell.broadinstitute.org/single_cell/data/public/SCP1830/spatial-atlas-of-molecular-cell-types-and-aav-accessibility-across-the-whole-mouse-brain?filename=metadata.csv 
#An account at the single-cell portal may be required.
#The metadata contains for each of their quantified cells the spatial location (i.e. which region/sub-region it was found in). We will aggregate that information into regions that are 
#equivalents of the regions that we analysed with our bulkseq data. We have performed the annotation of (sub-)structures in the Shi et al data and 'our' regions by hand based on the 
#method description in the Shi et al., manuscript.

#First we will load the metadata from Shi et al.
metadata_Wanglab <- read_csv("input_data/Shi_et al_STARmap_cellcounts/Shi_et_al_metadata.csv")
metadata_Wanglab <- metadata_Wanglab[-1,]
metadata_Wanglab <- as.data.frame(metadata_Wanglab)
#This dataset contains data from spinal cord and the brain. We will subset to the 'brain' data only. We will also focus here only on males, as data from across the brain is only reported
#well for males in the manuscript
table(metadata_Wanglab$organ__ontology_label)
metadata_Wanglab <- metadata_Wanglab %>% filter(organ__ontology_label == 'brain')
metadata_Wanglab <- metadata_Wanglab %>% filter(sex == 'male')
#Further, we will limit our analysis on the authors' sagittal data, to avoid potential disagreements between sagittal and coronal slices 
#We use the prefix in the NAME slot to pull out the sagittal data
metadata_Wanglab$dataset <- unlist(lapply(strsplit(as.character(metadata_Wanglab$NAME), split = "_"), function(x) x[[1]]))
metadata_Wanglab$cellID <- unlist(lapply(strsplit(as.character(metadata_Wanglab$NAME), split = "sagittal"), function(x) x[[2]]))

#We can obtain information about the number of 'tissues'/regions and the number of cells of each cell type found in these
sort(table(metadata_Wanglab$Tissue_Symbol), decreasing = T)
sort(table(metadata_Wanglab$Maintype_Symbol), decreasing = T)
#This appears to cover all relevant glia cell types, so we can continue working with this data
table(metadata_Wanglab$Maintype_Symbol, metadata_Wanglab$Tissue_Symbol)
#Now we load a 'dictionary' where we manually curated each of the sub-regions in the Shi et al dataset with regions that correspond to the ones we profiled via Bulk-Seq
Shi_Hahn_region_annotation <- read_delim("input_data/Shi_et al_STARmap_cellcounts/Shi_Hahn_region-annotation.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
#'region' represents the structures in Shi et al, while 'tissue_Ohahn' represents the manually-curated mapping of regions in our bulkseq data
#we use mapvalues to take the metadata's 'Tissue_Symbol' and relabel it using the informatino from the previously loaded annotation dataframe
metadata_Wanglab$CAregion <- plyr::mapvalues(metadata_Wanglab$Tissue_Symbol, 
                                             from = Shi_Hahn_region_annotation$region, 
                                             to = Shi_Hahn_region_annotation$tissue_Ohahn)

#We can inspect the results by using theh table function and check if cell types are now grouped into regions that have the same acronyms as our bulkseq data
table(metadata_Wanglab$Maintype_Symbol, metadata_Wanglab$CAregion)
#The outcome should look like this
#                                               cc   cer   cor    cp   ent    hi    hy   olf other   plx    th
# Astrocytes                                   1180  4956  6127  2136   344  4070  1761  3030 17367    56  2450
# Cerebellum neurons                              1 60393    10     1     0    12     2    10  2587     6     0
# Cholinergic and monoaminergic neurons          78   419   270   391     6   157   475   593  2572    20    60
# Choroid plexus epithelial cells                16    64    54    18     0    31     9    38   749  1804    20
# Dentate gyrus granule neurons                   0     2    22     0     3  4390     2     2    85     0     1
# Di- and mesencephalon excitatory neurons       43   314   209    87     9   105   340   321  6498     8  6778
# Di- and mesenphalon inhibitory neurons         87   125   113   106     1   125   235   381  5733     6    82
# Ependymal cells                                20   166    54    21     2    39    33   178  2459   110    11
# Glutamatergic neuroblasts                      86   461    85     8     5    34    13  9984   419     7    11

#We take those numbers to analyse the cell composition of each region/tissue
summary_stats <- as.data.frame.matrix((table(metadata_Wanglab$Maintype_Symbol, metadata_Wanglab$CAregion))) 

#We are interested in knowing what the composition of cells in a given region looks like
#For that purpose, we perform column-wise (region-wise) normalization
freqs_STARmap <- t(apply(t(summary_stats),1, function(x) x/sum(x)))
#We melt this into the long format. the resulting plot dataframe will be the basis for the following diagnostic plots
plot_df <- melt(freqs_STARmap)
colnames(plot_df) <- c("region", 'cell_type', 'abundance')
#There are several regions in the Shi dataset that we have no corresponding bulk tissue collected from (labelled 'other'). We exclude that
plot_df <- plot_df %>% filter(region != 'other')

#Now we want to compare the cell composition data with a region's CAS slope, that we quantified in the previous analysis
load('R_objects/BulkSeq_CAS_slopes.bin') #This loads the 'slope_dataframe' from the previous analysis

#Now we map over the regions' slope ('age.trend'). We do not have clear data in the Shi dataset for pons, medullar, svz. Further, we focus only one cortical (motor cortex) and hipppocampal are
#Visual cortex and posterior hippocampus were thus not inpsectec
plot_df$slope <-  as.numeric(as.character( plyr::mapvalues(plot_df$region, 
                                                           from = slope_dataframe$tissue,
                                                           to = slope_dataframe$age.trend)))
#Now we plot for each cell type its relative abundance across 10 regions, against the regions' CAS slopes
#This creates Figure S6
(myplot <- ggplot(plot_df, aes(x=slope, y=abundance, color=region)) + geom_point2()  + 
    scale_y_continuous(limits = c(0,NA), oob=squish) +
    scale_x_continuous(limits = c(0,0.02),expand = c(0,0)) + facet_wrap(~cell_type, scales = 'free_y') + geom_smooth(color='black', method = 'lm',  se = T)
)

#We will  also plot specifically the plots for the four major glia classes
#These create Figure 2J
(myplot <- ggplot(plot_df %>% filter(cell_type=='Oligodendrocyte precursor cells'), aes(x=slope, y=abundance, color=region)) + geom_point2()  + 
    scale_y_continuous(limits = c(0,0.04), oob=squish,expand = c(0,0)) +
    scale_x_continuous(limits = c(0.005,0.02),expand = c(0,0)) + facet_wrap(~cell_type, scales = 'free_y') + geom_smooth(color='black', method = 'lm',  se = T)
)

(myplot <- ggplot(plot_df %>% filter(cell_type=='Astrocytes'), aes(x=slope, y=abundance, color=region)) + geom_point2()  + 
    scale_y_continuous(limits = c(0,0.2), oob=squish,expand = c(0,0)) +
    scale_x_continuous(limits = c(0.005,0.02),expand = c(0,0)) + facet_wrap(~cell_type, scales = 'free_y') + geom_smooth(color='black', method = 'lm',  se = T)
)

(myplot <- ggplot(plot_df %>% filter(cell_type=='Microglia'), aes(x=slope, y=abundance, color=region)) + geom_point2()  + 
    scale_y_continuous(limits = c(0,0.06), oob=squish,expand = c(0,0)) +
    scale_x_continuous(limits = c(0.005,0.02),expand = c(0,0)) + facet_wrap(~cell_type, scales = 'free_y') + geom_smooth(color='black', method = 'lm',  se = T)
)

(myplot <- ggplot(plot_df %>% filter(cell_type=='Oligodendrocytes'), aes(x=slope, y=abundance, color=region)) + geom_point2()  + 
    scale_y_continuous(limits = c(0,0.7), oob=squish,expand = c(0,0)) +
    scale_x_continuous(limits = c(0.005,0.02),expand = c(0,0)) + facet_wrap(~cell_type, scales = 'free_y') + geom_smooth(color='black', method = 'lm',  se = T)
)

#Finally we will analyse for each cell type the data with spearman correlation, as well as a linear model.
#We will use the resulting pvalues to assess any potential CAS-cell type composition relationship
#First we iniate a dataframe in which we'll store the outputs
output_df <- data.frame( pvalue=1, cor_value=1, row.names = 'init')
for (cells in unique(plot_df$cell_type)) {
  #Perform correlation analysis
  results_df <- cor.test((plot_df %>% filter(cell_type==cells))$slope, (plot_df %>% filter(cell_type==cells))$abundance, method = 'spearman')
  output_df[cells, 'pvalue'] <- results_df$p.value
  output_df[cells, 'cor_value'] <- results_df$estimate
  #Build sub-dataset as input for linear regression
  plot_df_sub <- plot_df %>% filter(cell_type==cells)
  #perform linear model where CAS slope is explained by per-region cell type abundancd
  model <- lm(slope ~ abundance, data=plot_df_sub)
  model_summary <- summary(model)
  output_df[cells, 'Ftest_pvalue'] <- model_summary$coefficients[8]
  output_df[cells, 'lm_estimate'] <- model_summary$coefficients[2]
  
}
#Remove initial row
output_df <- output_df[-1, ]
#Perform multiple testing correction
output_df$padj <- p.adjust(output_df$pvalue, method = 'BH')
output_df$padj_lm <- p.adjust(output_df$Ftest_pvalue, method = 'BH')

#We find no significant association for any cell type. This provides evidence that there is no direct relationship between cell type composition and CAS change.
