###Here, we perorm the DEG analsyis for the YMP/PBS injection group
#First we load the Deseq2 object
load(file='R_objects/dds_BulkSeq_YMP.bin')

design(dds_CA2_YMP) <- ~ treatment + tissue
#Run Deseq2
dds_CA2_YMP <- DESeq(dds_CA2_YMP, fitType = 'local')
dds_CA2_YMP_list <- list(All=dds_CA2_YMP
)


#Now we create a results list onto which we will reiterativley add results tables from the Deseq2 analysis
results_list_CA2_YMP <- list()
results_list_CA2_YMP$Female <- list()

#Define cutoff for differntial expression test 
padj_cutoff <- 0.05
#Set the 'sex' that hsould be analysed in this run - we only have male YMP/PBS mice, so we stick with that
sex <- 'Male'

#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective region and then run the regular Deseq2 workflow
for (tissue in c( unique(as.character(dds_CA2_YMP_list$All$tissue)))) {
  print(tissue)
  #subset for current tissue and re-run Deseq2
  if (tissue != 'All') {
    dds_CA2_YMP_list[[tissue]] <- dds_CA2_YMP[, dds_CA2_YMP$tissue == tissue]
    dds_CA2_YMP_list[[tissue]]$tissue <- droplevels(dds_CA2_YMP_list[[tissue]]$tissue)
    design(dds_CA2_YMP_list[[tissue]]) <- ~treatment
    #Run Deseq2 and store in the Deseq2 object list
    dds_CA2_YMP_list[[tissue]] <- DESeq(dds_CA2_YMP_list[[tissue]], fitType = 'local')
  }
  
  #Now take the Deseq2 object and run the extraction of the pairwise comparisons
  dds_tissue_temp <- dds_CA2_YMP_list[[tissue]]
  results_list_tissue <- list()
  
  #Create a list of all pairwise comparisons
  comparison_list <- t(combn(unique(colData(dds_tissue_temp)$treatment),2))
  
  #Iterate over all pairwise comparisons and extract the resultstable, store it in the results list. then flip the conditions (so that results for 3 vs 21 and 21 vs 3 months is stored)
  for (row in 1:nrow(comparison_list)) {
    print(row)
    #Get the conditions/timepoints to test
    cond1 <- as.character(comparison_list[row,1])
    cond2 <- as.character(comparison_list[row,2])
    
    folder <- paste(cond1, cond2, sep = "_vs_")
    print(folder)
    
    aspect_to_test <- 'treatment'
    
    results_temp <- results(dds_tissue_temp, contrast=c(aspect_to_test, cond1,cond2), cooksCutoff=T) #Set contrast and retrieve results
    results_temp$gene_symbol <- row.names(results_temp) #get the rownames as a separate column in the results report
    resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),]) #create a simple dataframe ordered by the padj
    resSig <- subset(resOrdered, padj < padj_cutoff) #create a simple dataframe containing only the significant results
    dim(resSig)
    #Store the output tables in the results_list
    results_list_tissue[[folder]]$resall <- results_temp
    results_list_tissue[[folder]]$resOrdered <- resOrdered
    results_list_tissue[[folder]]$ressig <- resSig
    
    
    ##Flip conditions and re-run the results extraction
    
    cond1 <- as.character(comparison_list[row,2])
    cond2 <- as.character(comparison_list[row,1])
    folder <- paste(cond1, cond2, sep = "_vs_")
    print(folder)
    
    dds_tissue_temp <- dds_tissue_temp
    
    aspect_to_test <- 'treatment'
    
    results_temp <- results(dds_tissue_temp, contrast=c(aspect_to_test, cond1,cond2), cooksCutoff=T) #Set contrast and retrieve results
    results_temp$gene_symbol <- row.names(results_temp) #get the rownames as a separate column in the results report
    resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),]) #create a simple dataframe ordered by the padj
    resSig <- subset(resOrdered, padj < padj_cutoff) #create a simple dataframe containing only the significant results
    dim(resSig)
    #Store the output tables in the results_list
    results_list_tissue[[folder]]$resall <- results_temp
    results_list_tissue[[folder]]$resOrdered <- resOrdered
    results_list_tissue[[folder]]$ressig <- resSig
    
    
    
  }
  #Store tissue results in the major results list, then save the ouput and start the loop for the next tissue
  results_list_CA2_YMP[[sex]][[tissue]]  <- results_list_tissue
  #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
  save(dds_CA2_YMP_list, results_list_CA2_YMP,  file='R_objects/dds_BulkSeq_Rejuvenation_YMP.bin')
  
}


#Having completed the analysis, we're intersted in checking how many DEGs we found in each region, 
#so we create a bargraph of DEGs for each region akin to the graphs for the aging cohort in Figure S4A


#We will also create a bargraph that summarises all DEGs per regions (passing the padj cutoff in at least two comparisons with 3m)
#We set a padj cutoff of 0.05 to mark a gene as differentially expressed.
padj_cutoff <- 0.05

plot_df_crossTissue <- data.frame(direction='A', Comparison='A', Numb_of_genes=1, tissue='A', sex='A')
#We will focus here on the data that covers both sexes
results_list_CA2_YMP_Male <- results_list_CA2_YMP$Male
#We'll iterate over each tissue, collapse the results lists from all the age-related comparisons and count in how many of these 
#a given gene passed the significance threshold

for (tissue in names(results_list_CA2_YMP_Male)) {
  
  print(paste('Processing tissue:', tissue, sep = ' '))
  #We subset for the results tables from the current tissue
  results_list_tissue <- results_list_CA2_YMP_Male[[tissue]]
  #We create a dataframe that will store the information 
  barplot_frame <- data.frame()
  
  #We will need to know in which direction a given gene is moving (up/downregualted)
  #To this end, we'll take the YMP vs PVS-comparison as reference and use the sign of the log2FC as indicator of the regulation direction
  regulation_df <- results_list_tissue[["YMP_vs_PBS"]]$resOrdered
  #We use lappy to extract the DEGs
  DEGs_of_tissue <- (unlist((lapply(results_list_tissue[c('YMP_vs_PBS' )], function (x) (row.names(subset(x[["resall"]], padj < padj_cutoff)))))))
  #now we count often a given gene 'occured'  in the DEGS_of_tissue vector and supset for genes that have been found at least twice (i.e. passing two times the significance threshould)
  DEGs_of_tissue <- names(which(table(DEGs_of_tissue) >= 1))
  #We'll take the DEGs and subset the regulation_df for these
  regulation_df <- regulation_df[DEGs_of_tissue,]
  #From the log2FC in that data frame we extract the number of up/downregulated genes
  down_regulated <- nrow(regulation_df %>% filter(log2FoldChange < 0))
  up_regulated <- nrow(regulation_df  %>% filter(log2FoldChange > 0))
  #and now we'll store that information in the barplot frame 
  barplot_frame["down","TotalDEG"] <- data.frame(down_regulated)
  barplot_frame["up","TotalDEG"] <- data.frame(up_regulated)
  #Now we'll transform that to make it easier for ggplot to work with it down the road
  plot_df <- melt(as.matrix(barplot_frame))
  colnames(plot_df) <- c('direction',"Comparison", "Numb_of_genes")
  #We also provide the information about the current tissue
  plot_df$tissue <- tissue
  plot_df$sex <- 'Both'
  plot_df_crossTissue <- rbind(plot_df_crossTissue, plot_df)
  
  
}

plot_df_crossTissue <- plot_df_crossTissue[-1,] #Remove the first row we needed to set up the dataframe
#we'll set the order of factors in the regulation column 
plot_df_crossTissue$direction <- factor(plot_df_crossTissue$direction, levels = c("up","down"), ordered = T)
#We'll turn the tissue/region column into a vector with factors ordered according to the total number of DEGs
plot_df_crossTissue$tissue <- factor(plot_df_crossTissue$tissue, levels = rev(as.character((plot_df_crossTissue %>% dplyr::group_by(tissue)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
                                     ordered = T)

#Now we plot this and generate Figure 7F
ggplot(plot_df_crossTissue , aes(x=tissue, y=Numb_of_genes, fill=direction), color=NA) + geom_bar(stat = "identity")


