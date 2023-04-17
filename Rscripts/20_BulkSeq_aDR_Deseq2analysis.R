###Here, we perorm the DEG analsyis for the aDR/AL treatment group
#First we load the Deseq2 object
load(file='R_objects/dds_BulkSeq_DR.bin')
#Dietary restriction induces a known differential exprssion in the liver, which we can use to confirm that the DR protocol was effective, so we will focus on the liver dataset first
dds_CA2_DR_liver  <- dds_CA2_DR[, dds_CA2_DR$tissue == "LL"]
dds_CA2_DR_liver$treatment <- droplevels(dds_CA2_DR_liver$treatment)
dds_CA2_DR_liver$tissue <- droplevels(dds_CA2_DR_liver$tissue)
design(dds_CA2_DR_liver) <- ~ treatment
#Run Deseq2
dds_CA2_DR_liver <- DESeq(dds_CA2_DR_liver, fitType = 'local')

#We'll look at the overall effect of DR on the transcriptome via PCA. We will use the varianceStabilizingTransformation on the count matrix and focus on the top 5000 expressed genes
#This creates panel S15C
plotPCA(varianceStabilizingTransformation(dds_CA2_DR_liver), intgroup = "treatment", ntop=5000)
results_table_DR_Liver <- as.data.frame(results(dds_CA2_DR_liver, contrast=c('treatment', 'DR', 'AL'), cooksCutoff=T))
results_table_DR_Liver$gene_symbol <- row.names(results_table_DR_Liver)
resSig <- subset(results_table_DR_Liver, padj < 0.05)

#We'll check the number of differentially regulated genes
#This creates panel S15D
table(resSig$log2FoldChange > 0)

#This looks good. We see that aDR leads to significant expression changes and seperates the diet groups in the PCA.
#We can also compare these changes to published DEG results from liver tissue of 3 months-long aDR-fed mice in the B6D2F1 background
#We have previously published these results and can load&comapre these

B6D2F1_aDR_vs_AL_liverDeseq2results <- read_delim("input_data/BulkSeq_Rejuvenation/Hahn_et_al_2019_aDR_vs_AL_liverDeseq2results.txt", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
B6D2F1_aDR_vs_AL_liverDeseq2results <- as.data.frame(B6D2F1_aDR_vs_AL_liverDeseq2results)

#We'll check seprately for up- and down-regualted genes if and to what degree DEGs overlap
#First, we'll check the number of genes being up-regulated in both datasets
#We extract the DEGs from each dataset individually that have a postitive (up-regulated under aDR) log2FoldChange
gene_list_in_HahnEtAl <- B6D2F1_aDR_vs_AL_liverDeseq2results %>% filter(log2FoldChange > 0) %>% filter(padj <= 0.05) %>% pull(gene_symbol)

gene_list_in_LiverOfDRmice <- results_table_DR_Liver %>% filter(log2FoldChange > 0)  %>% filter(padj <= 0.05) %>% pull(gene_symbol)
#Then we do run intersect to get the common DEGs
common_genes <- intersect(gene_list_in_HahnEtAl,
                          gene_list_in_LiverOfDRmice)
#And then we need a list of expressed genes, for which we use the current liver dataset. We use Deseq2's metric for independent filtering, which is marked by genes 
#exhibiting no 'na' in the padj column
Expressed_in_tissue <- results_table_DR_Liver %>% filter(!is.na(padj)) %>% pull(gene_symbol)
#now we build a 2x2 matrix out of these gene lists as input for Fisher's test
m <- matrix(c(length(common_genes), length(gene_list_in_HahnEtAl) - length(common_genes),  length(gene_list_in_LiverOfDRmice) - length(common_genes), length(Expressed_in_tissue) - length(gene_list_in_LiverOfDRmice) - length(gene_list_in_HahnEtAl) - length(common_genes)), ncol = 2)
#And run this 
fisher.test(m)
#A highly significant outcome. The genes are regulated to a significant degree in the same direciton
#We can use vennerable's Venn function to create Figure S15E
Vestem <- Venn(list(Hahn=gene_list_in_HahnEtAl, 
                    Liver=gene_list_in_LiverOfDRmice)
)

plot(Vestem, doWeights = T)
#This analysis can be repeated for the down-regualted genes to obtain the second half of the Figure S15E


#Having confirmed the validity of our aDR protcol, we can now analyse the bulkseq data from the brain
load(file='R_objects/dds_BulkSeq_DR.bin')
#First, we'll filter out the liver samples and cleanup the factor levels
dds_CA2_DR  <- dds_CA2_DR[, !dds_CA2_DR$tissue == "LL"]
dds_CA2_DR$treatment <- droplevels(dds_CA2_DR$treatment)
dds_CA2_DR$tissue <- droplevels(dds_CA2_DR$tissue)
design(dds_CA2_DR) <- ~ treatment + tissue
#Run Deseq2
dds_CA2_DR <- DESeq(dds_CA2_DR, fitType = 'local')
dds_CA2_DR_list <- list(All=dds_CA2_DR
)


#Now we create a results list onto which we will reiterativley add results tables from the Deseq2 analysis
results_list_CA2_DR <- list()
results_list_CA2_DR$Female <- list()

#Define cutoff for differntial expression test 
padj_cutoff <- 0.05
#Set the 'sex' that hsould be analysed in this run - we only have female aDR/AL mice, so we stick with that
sex <- 'Female'

#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective region and then run the regular Deseq2 workflow
for (tissue in c( unique(as.character(dds_CA2_DR_list$All$tissue)))) {
  print(tissue)
  #subset for current tissue and re-run Deseq2
  if (tissue != 'All') {
    dds_CA2_DR_list[[tissue]] <- dds_CA2_DR[, dds_CA2_DR$tissue == tissue]
    dds_CA2_DR_list[[tissue]]$tissue <- droplevels(dds_CA2_DR_list[[tissue]]$tissue)
    design(dds_CA2_DR_list[[tissue]]) <- ~treatment
    #Run Deseq2 and store in the Deseq2 object list
    dds_CA2_DR_list[[tissue]] <- DESeq(dds_CA2_DR_list[[tissue]], fitType = 'local')
  }
  
  #Now take the Deseq2 object and run the extraction of the pairwise comparisons
  dds_tissue_temp <- dds_CA2_DR_list[[tissue]]
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
  results_list_CA2_DR[[sex]][[tissue]]  <- results_list_tissue
  #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
  save(dds_CA2_DR_list, results_list_CA2_DR,  file='R_objects/dds_BulkSeq_Rejuvenation_DR.bin')
  
}


#Having completed the analysis, we're intersted in checking how many DEGs we found in each region, 
#so we create a bargraph of DEGs for each region akin to the graphs for the aging cohort in Figure S4A


#We will also create a bargraph that summarises all DEGs per regions (passing the padj cutoff in at least two comparisons with 3m)
#We set a padj cutoff of 0.05 to mark a gene as differentially expressed.
padj_cutoff <- 0.05

plot_df_crossTissue <- data.frame(direction='A', Comparison='A', Numb_of_genes=1, tissue='A', sex='A')
#We will focus here on the data that covers both sexes
results_list_CA2_DR_Female <- results_list_CA2_DR$Female
#We'll iterate over each tissue, collapse the results lists from all the age-related comparisons and count in how many of these 
#a given gene passed the significance threshold

for (tissue in names(results_list_CA2_DR_Female)) {
  
  print(paste('Processing tissue:', tissue, sep = ' '))
  #We subset for the results tables from the current tissue
  results_list_tissue <- results_list_CA2_DR_Female[[tissue]]
  #We create a dataframe that will store the information 
  barplot_frame <- data.frame()
  
  #We will need to know in which direction a given gene is moving (up/downregualted)
  #To this end, we'll take the aDR vs AL-comparison as reference and use the sign of the log2FC as indicator of the regulation direction
  regulation_df <- results_list_tissue[["aDR_vs_AL"]]$resOrdered
  #We use lappy to extract the DEGs
  DEGs_of_tissue <- (unlist((lapply(results_list_tissue[c('aDR_vs_AL' )], function (x) (row.names(subset(x[["resall"]], padj < padj_cutoff)))))))
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

#Now we plot this and generate Figure 7B
ggplot(plot_df_crossTissue , aes(x=tissue, y=Numb_of_genes, fill=direction), color=NA) + geom_bar(stat = "identity")


