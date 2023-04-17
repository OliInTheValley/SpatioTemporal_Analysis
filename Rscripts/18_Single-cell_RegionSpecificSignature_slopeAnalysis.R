
##Here we will plot and quantify the trajectories of region-specific aging signatures iteratively cell types from a single-cell/nuc dataset, in this case from the aging nuc-seq hippocampus dataset
#We are trying to investigate which cell type is predominantly associated with the increase of the region-specific signatures, i.e. the aging trajectory of which cell type is mainly represented by the signature found in bulk
#To do that, we will iterate through each region-specific signature, retrieve the per-cell type scores and perform the differential slope analysis akin to what we did previously with the CAS. 
#That is, we will build a model which we can use to test slopes from multiple tissues in comparison
#We will collect that information for each signature and visualize the output as a heatmap
#We deomnstrate this analysis here using the example of the hippocampus. but the same analysis can be apllied to the caudate putamen dataset
#First we load the output of the VISION score calcualtion from the nuc-seq analysis
load('R_objects/Nucseq_Hip_AgingDEGs_Scoretable.bin')


####We will make a copy of the original score table. 
score_table <- scRNA_Aging_Scoretable 
#We make sure 'age' is a numerical vector
score_table$age <- plyr::mapvalues(score_table$age, from =  c('Young','Old'), to = c(3,21))
score_table$age <-  as.numeric(as.character(score_table$age))

#We'll assign the rownames into a new column called 'sampleID', which we'll need in the loop below

#We have chosen a couple of region-specific signatures that we previously confirmed to be 'truely' region-specific, and that represent a range of region bioloigy: cortical, hippocampal, cerebellar, striatal and 'white matter' 
region_signatures <- c('cp', 'hi', 'cc', 'cer', 'cor', 'vis')
cell_types_to_analyze <- as.character(unique(score_table$cell_type))

#We will  run models across all regions together, similar to what we did with the CAS before - to perform differential slope analysis. We will need respective dataframes to hold that informatino too
#We note that we forgo a seperate analysis where we would assess specifically if a signature in a given cell type would associate with age, as the high number of cells/degrees of freedom render the anlaysis useless, as essentially eveything is significant
slope_dataframe_collectionDF <- data.frame(init=rep(1, length(cell_types_to_analyze)), row.names = cell_types_to_analyze)
pvalue_dataframe_collectionDF <- data.frame(init=rep(1, length(cell_types_to_analyze)), row.names = cell_types_to_analyze)
stars_dataframe_collectionDF <- data.frame(init=rep('*', length(cell_types_to_analyze)), row.names = cell_types_to_analyze)



#We also create some results lists to store model coefficients and results from the differential slope analysis
cofficent_frame_list <- list() 
slope_PairWiseResults_list <- list() 

#Now we start our loop - it will run for each region-specific signature
for (tissue_score in region_signatures) {
  print(tissue_score)
  
  #To keep it simple, we'll just assign this column the label 'score'
  score_table$score <- score_table[,paste('Unique_', tissue_score,'_UpDownReg_AgingSignature', sep = '')]
  
  #Now, we will calculate a linear model that allows for an interaction between age and cell type (denoted with an)
  m.interaction <- lm(score ~ age*cell_type, data = score_table)
  #We are going to setup a dataframe that contains the intercepts for each cell type for later inspection of the offset
  coefficient_dataframe <- as.data.frame(m.interaction$coefficients)
  colnames(coefficient_dataframe) <- 'coefficients'
  coefficient_dataframe <- coefficient_dataframe[grep('^cell_type', row.names(coefficient_dataframe), value = T),,drop=F]
  coefficient_dataframe$cell_type <- gsub('cell_type', '',row.names(coefficient_dataframe))
  colnames(coefficient_dataframe)[1] <- 'offset'
  
  #Now we'll use the function lstrends (fromt he lsmeans package) to quantify slopes with confidence intervals to allow cell type-to-cell type comparison
  m.lst <- lstrends(m.interaction, "cell_type", var="age")
  #Let's look at that output and store it in a separte dataframe to hold the information about the slopse
  slope_dataframe <- as.data.frame(m.lst)
  
  #Create pairwise comparisons between slopes and test using the pais function and store that result
  #This lets us know which slopes are statistically significantly different from one anothre
  slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))
  
  #We will now plot the results from lstrends as bargraphs with confidence intervals
  #We'll order the rank of the cell types in ascending order along their cell type-specific score slope estimate
  slope_dataframe$cell_type <- factor(slope_dataframe$cell_type, levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = T),'cell_type']),
                                   ordered = T)
  #We will color the bars accoring to the score incrase ('age.trend')
  #If the region-signature is in fact representing a specific cell type's aging trajectory, we should see that specific cell type(s) stand outperform other cell types substantially
  (myplot <- ggplot(slope_dataframe, aes(x=cell_type, y=age.trend,fill=age.trend)) +  
      geom_bar(position=position_dodge(), stat="identity", mapping = aes(fill=age.trend)) + xlab("") + ylab("CAS slope")  +
      scale_fill_gradient2(low=brewer.pal(11, "RdYlBu")[10], midpoint = quantile(slope_dataframe$age.trend)[3], mid=brewer.pal(11, "RdYlBu")[6], high = brewer.pal(11, "RdYlBu")[2])+
      geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      scale_y_continuous(expand = c(0,0)) + 
      theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1))
  )
  
  #Prepare coeifficent_dataframe to assemble with slopes for this score with slopes from other scores
  slope_dataframe$signature <- tissue_score
  row.names(slope_dataframe) <- slope_dataframe$cell_type
  cofficent_frame_list[[tissue_score]] <- slope_dataframe 
  slope_PairWiseResults_list[[tissue_score]] <- slope_pair_wiseTestResults
  #Now we're mainly interested in the slope of the cell type with the highest slope (as this could be the 'candidate') So we define that cell type 
  #of interest as 'reference' and grab only the pariwise comparisons to that particular region 
  reference_cell <- as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = T),'cell_type'])[1]
  slope_pair_wiseTestResults_onlyRef <- slope_pair_wiseTestResults[grepl(reference_cell, slope_pair_wiseTestResults$contrast),]
  
  #We'll reformat this so that we can plot this information on a heatmap
  p_valueDf <- data.frame(contrast='reference', estimate=1, SE=1, df=1, t.ratio=1, p.value=2, row.names = reference_cell)
  for (current_celltype in everything_but(cell_types_to_analyze, c(reference_cell)) ){ 
    new_line <- slope_pair_wiseTestResults_onlyRef[grepl(current_celltype, slope_pair_wiseTestResults_onlyRef$contrast),]
    row.names(new_line) <- current_celltype
    p_valueDf <- rbind(p_valueDf,new_line)
  }
  #and create a version of it with astericks
  p_valueDf$stars <- apply(p_valueDf['p.value'], 1, stars.return)
  p_valueDf[reference_cell, 'stars'] <- 'ref'
  
  #We store the information in the respective dataframe
  slope_dataframe_collectionDF[cell_types_to_analyze, tissue_score] <-   slope_dataframe[cell_types_to_analyze,'age.trend']
  pvalue_dataframe_collectionDF[cell_types_to_analyze, tissue_score] <-   p_valueDf[cell_types_to_analyze,'p.value']
  stars_dataframe_collectionDF[cell_types_to_analyze, tissue_score] <-   p_valueDf[cell_types_to_analyze,'stars']
  
  
}
#We reomve the initial column that we needed when setting up the dataframe
slope_dataframe_collectionDF$init <- NULL
stars_dataframe_collectionDF$init <- NULL

#Now we're putting it all together
#We want to plot the aggreagted informatino as a heatmap, where columns are the quantified per-cell type slopes, and rows are the signatures
#Each cell thus represents the slope for a given signature in a given cell type
#first we remove the first column which we needed to set up the dataframe in the first place
pvalueIndividual_dataframe_collectionDF$init <- NULL

colfunc <- colorRampPalette(c("#19F5FF","white", "red"))

max_value <- 1

#This creates panels 6H/K
pheatmap(t(slope_dataframe_collectionDF), 
         show_colnames = T, 
         cluster_cols = T, cluster_rows = T,  
         clustering_method = "single", border_color = NA,display_numbers = t(stars_dataframe_collectionDF) )

