
##Here we will plot and quantify the trajectories of region-specific aging signatures iteratively across each region
#We aim to assess if region-specific signatures really reflect region-specific dynamics. 
#To do that, we will iterate through each region-specific signature check in which regions the score associates with age and
#then do the differential slope analysis akin to what we performed previously with the CAS. That is, we will build a model which we can use to test slopes from multiple tissues in comparison
#We will collect that information for each signature and visualize the output as a heatmap
#First we load the output of the VISION score calcualtion
load('R_objects/BulkSeq_Aging_Scoretable.bin')


####We will make a copy of the original score table. 
score_table <- BulkSeq_Aging_Scoretable 
#We'll assign the rownames into a new column called 'sampleID', which we'll need in the loop below
score_table$sampleID <- row.names(score_table)

#We did not find a convincing set of genes to form a signature for the entorhinal cortex, so we'll filter that tissue out
#We also did not build a signautre for the posterior hippocampus, as we did not want to have the two very similiar hippocampal regions 'compete' for the dectection of regon-specific genes
#So we'll filter that one out too
score_table <- score_table %>% filter(tissue != 'hi2') %>% filter(tissue != 'ent')

#Now we retrieve the the names of regions for which we found and quantified a signature - we could do this by hand but this is a better route to avoid calling a column that doesnt exist
scores_to_plot <- grep('Unique_[a-z]*_UpDownReg_AgingSignature',colnames(score_table), value = T)
tissues_with_scores <- unlist(lapply(strsplit(scores_to_plot, split = "_"), function(x) x[[2]]))

#We create dataframes that will store information about slope, pvalue and pvalue-converted-to-asterisks for later plotting
#In this analysis we will run models on each region individually - to test if they associate with age at all. So we need dataframes to hold the information from those outputs
slopeIndividual_dataframe_collectionDF <- data.frame(init=rep(1, length(tissues_with_scores)), row.names = tissues_with_scores)
pvalueIndividual_dataframe_collectionDF <- data.frame(init=rep(1, length(tissues_with_scores)), row.names = tissues_with_scores)
starsIndividual_dataframe_collectionDF <- data.frame(init=rep('*', length(tissues_with_scores)), row.names = tissues_with_scores)

#And we will also run models across all regions together, similar to what we did with the CAS before - to perform differential slope analysis. We will need respective dataframes to hold that informatino too
slope_dataframe_collectionDF <- data.frame(init=rep(1, length(tissues_with_scores)), row.names = tissues_with_scores)
pvalue_dataframe_collectionDF <- data.frame(init=rep(1, length(tissues_with_scores)), row.names = tissues_with_scores)
stars_dataframe_collectionDF <- data.frame(init=rep('*', length(tissues_with_scores)), row.names = tissues_with_scores)


#We also create some results lists to store model coefficients and results from the differential slope analysis
cofficent_frame_list <- list() 
slope_PairWiseResults_list <- list() 

#Now we start our loop - it will run for each region-specific signature
for (tissue_score in tissues_with_scores) {
  print(tissue_score)
  
  #To keep it simple, we'll just assign this column the label 'score'
  score_table$score <- score_table[,paste('Unique_', tissue_score,'_UpDownReg_AgingSignature', sep = '')]
  #We make sure 'age' is a numerical vector
  score_table$age <-  as.numeric(as.character(score_table$age))
  
  #First, we'll assess in each region inddividually, if the region is associated with age. We call those the 'individual models'
  #Run a per tissue linear regression for comparisons
  #We build a temporary datafrate that we need to initialize
  individual_slope_dataframe <- data.frame(tissue='init', age.trend=0, p.value=1, row.names = 'init', stringsAsFactors = F)
  for (tissue_to_analyze in tissues_with_scores) {
    #We'll buid a linear model based on data that is subsetted to a given region
    m.tissue <- lm(score ~ age, data = score_table %>% filter(tissue == tissue_to_analyze))
    results_df <- summary(m.tissue)
    #Then we extract the relevant information regarding the slope and the p value for the F-statistic
    individual_slope_dataframe[tissue_to_analyze, 'tissue'] <- tissue_to_analyze
    individual_slope_dataframe[tissue_to_analyze, 'age.trend'] <- results_df$coefficients[2]
    individual_slope_dataframe[tissue_to_analyze, 'p.value'] <- results_df$coefficients[8]
  }
  individual_slope_dataframe <- individual_slope_dataframe[-1,]
  #We corrct the pvalues for multiple testing
  individual_slope_dataframe$p.adj <- p.adjust(individual_slope_dataframe$p.value, method = 'bonferroni')
  #And then convert those into astericks to make it easier to inspect the results
  individual_slope_dataframe$stars <- apply(individual_slope_dataframe['p.adj'], 1, stars.return)

  #Having completed that, we'll store the results for this signature   
  slopeIndividual_dataframe_collectionDF[tissues_with_scores, tissue_score] <-   individual_slope_dataframe[tissues_with_scores,'age.trend']
  pvalueIndividual_dataframe_collectionDF[tissues_with_scores, tissue_score] <-   individual_slope_dataframe[tissues_with_scores,'p.adj']
  starsIndividual_dataframe_collectionDF[tissues_with_scores, tissue_score] <-   individual_slope_dataframe[tissues_with_scores,'stars']
  
  ###Now we move over to the cross-region analysis. this is aking to what we did previously for the CAS
  #We acknowledge that we could just drop all regions that do not show a significant association with age in the above's model but we'd rather filter the results later, as we'd create a different
  #statistical basis for each signature. E.g. we'd have differenig tissues we'd have to perform multiple testing on etc. So we run the cross-slope analysis across all regions and filter out the final results
  
  #First we plot the scores for each region along the age-axis. We also quanitfy the per-age group in a boxplot
  #This can create panels like the on in figure 6B for the caudate putmaen signature
  (myplot <- ggplot(score_table, aes(x=age, y=score)) + geom_boxplot(notch = F, outlier.colour = NA,  mapping = aes(fill=as.factor(age))) + 
      geom_smooth(data =  score_table, mapping = aes(x=age, y=score), color='red', se = T,method = 'lm') +
      geom_point(color='#bbbdbf', position = 'jitter') + facet_grid(~tissue) 
  )
  
  ####Slope analysis
  
  #Now we just focus on the trajectories, leaving out the faceting, boxplots and individual points
  #Because loess models are challengnign to interpret, we'll recreate the same plot using linear models instead
  (myplot <- ggplot(score_table, aes(x=age, y=score, color=tissue)) + 
      geom_smooth(data = score_table, mapping = aes(x=age, y=score), se = T,  method = 'lm') + scale_y_continuous(expand = c(0,0))
  )

  
  #Now, we will calculate a linear model that allows for an interaction between age and tissue (denoted with an)
  m.interaction <- lm(score ~ age*tissue, data = score_table)
  #We are going to setup a dataframe that contains the intercepts for each tissue for later inspection of the offset
  coefficient_dataframe <- as.data.frame(m.interaction$coefficients)
  colnames(coefficient_dataframe) <- 'coefficients'
  coefficient_dataframe <- coefficient_dataframe[grep('^tissue', row.names(coefficient_dataframe), value = T),,drop=F]
  coefficient_dataframe$tissue <- gsub('tissue', '',row.names(coefficient_dataframe))
  colnames(coefficient_dataframe)[1] <- 'offset'
  
  #Now we'll use the function lstrends (fromt he lsmeans package) to quantify slopes with confidence intervals to allow tissue-to-tissue comparison
  m.lst <- lstrends(m.interaction, "tissue", var="age")
  #Let's look at that output and store it in a separte dataframe to hold the information about the slopse
  slope_dataframe <- as.data.frame(m.lst)
  #   tissue    age.trend           SE  df      lower.CL    upper.CL
  #1      cp 0.0032433121 0.0005008826 712  2.259929e-03 0.004226696
  #2     cer 0.0007644110 0.0004885839 712 -1.948265e-04 0.001723648
  #3     plx 0.0034424847 0.0004897150 712  2.481027e-03 0.004403943
  #4      cc 0.0046792611 0.0004972849 712  3.702941e-03 0.005655581
  #Create pairwise comparisons between slopes and test using the pais function and store that result
  #This lets us know which slopes are statistically significantly different from one anothre
  slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))
  
  #We will now plot the results from lstrends as bargraphs with confidence intervals
  #We'll order the rank of the tissues in ascending order along their tissue-specific score slope estimate
  slope_dataframe$tissue <- factor(slope_dataframe$tissue, levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = T),'tissue']),
                                   ordered = T)
  #We will color the bars accoring to the score incrase ('age.trend')
  #This will create panels like the one shown in Figure 6C and Supplementary Figures 12 and 13
  #Ideally, we should see that the region where we learned the region-specific signature from, should be the one with the highest slope and in the best case, outperform each other region substantially
  (myplot <- ggplot(slope_dataframe, aes(x=tissue, y=age.trend,fill=age.trend)) +  
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
  row.names(slope_dataframe) <- slope_dataframe$tissue
  cofficent_frame_list[[tissue_score]] <- slope_dataframe 
  slope_PairWiseResults_list[[tissue_score]] <- slope_pair_wiseTestResults
  #Now we're mainly interested in the slope of the region where we learned the signature from (let's say caudate putamen) performs in comparison to all other regions. So we define that region 
  #of interest as 'reference' and grab only the pariwise comparisons to that particular region 
  slope_pair_wiseTestResults_onlyRef <- slope_pair_wiseTestResults[grepl(tissue_score, slope_pair_wiseTestResults$contrast),]
  
  #We'll reformat this so that we can plot this information on a heatmap
  p_valueDf <- data.frame(contrast='reference', estimate=1, SE=1, df=1, t.ratio=1, p.value=2, row.names = tissue_score)
  for (tissue in everything_but(tissues_with_scores, c(tissue_score)) ){ 
    new_line <- slope_pair_wiseTestResults_onlyRef[grepl(tissue, slope_pair_wiseTestResults_onlyRef$contrast),]
    row.names(new_line) <- tissue
    p_valueDf <- rbind(p_valueDf,new_line)
  }
  #and create a version of it with astericks
  p_valueDf$stars <- apply(p_valueDf['p.value'], 1, stars.return)
  p_valueDf[tissue_score, 'stars'] <- 'ref'
  
  #We store the information in the respective dataframe
  slope_dataframe_collectionDF[tissues_with_scores, tissue_score] <-   slope_dataframe[tissues_with_scores,'age.trend']
  pvalue_dataframe_collectionDF[tissues_with_scores, tissue_score] <-   p_valueDf[tissues_with_scores,'p.value']
  stars_dataframe_collectionDF[tissues_with_scores, tissue_score] <-   p_valueDf[tissues_with_scores,'stars']
  

}
#We reomve the initial column that we needed when setting up the dataframe
slope_dataframe_collectionDF$init <- NULL
stars_dataframe_collectionDF$init <- NULL

#Now we're putting it all together
#We want to plot the aggreagted informatino as a heatmap, where columns are the quantified slopes, and rows are the signatures
#Each cell thus represents the slope for a given signature in a given tissue
#However, we first want to eliminate each cell where the region's score did not associate with age in the first place - which we identified in the first part of the loop
#
#For that we take the pvalue dataframe where we collected from the individual tissue-score models
#first we remove the first column which we needed to set up the dataframe in the first place
pvalueIndividual_dataframe_collectionDF$init <- NULL
#Further, since we stored the information colum-wise, we'll need to transpose the dataframe
pvalueIndividual_dataframe_collectionDF <- as.data.frame(t(pvalueIndividual_dataframe_collectionDF))

#Then we filter out in the slope dataframe every cell where the individual models found no significant association with age 
slope_dataframe_collectionDF[pvalueIndividual_dataframe_collectionDF >=  0.05] <- NA 
#Same thing for the stats for the differential slope comparisons
#But first, we'll also have to transpose that one
stars_dataframe_collectionDF <- as.data.frame(t(stars_dataframe_collectionDF))
stars_dataframe_collectionDF[pvalueIndividual_dataframe_collectionDF >=   0.05] <- '' 
stars_dataframe_collectionDF <- as.data.frame(stars_dataframe_collectionDF)
colfunc <- colorRampPalette(c("#19F5FF","white", "red"))

max_value <- 1

#This creates panel S13 B 
pheatmap(slope_dataframe_collectionDF, 
         show_colnames = T, 
         cluster_cols = F, cluster_rows = F,  
         clustering_method = "ward.D2", border_color = NA, display_numbers = stars_dataframe_collectionDF)

