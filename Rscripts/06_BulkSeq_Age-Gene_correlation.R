
rm(list=ls())

##Here, we perform will correlate expression for each gene in each tissue/region with age and quantify the number of significant age-correlated genes


#First, we load our Deseq2 object as it contains all the counts and metadata
load('R_objects/dds_BulkSeq_Aging.bin')


#We'll setup a results list in which we'll store the output of the following analysis
BulkSeq_Agecorrelation_results_list <- list()

#We will rotate over each tissue/region (except the 'All') object
for (tissue in everything_but(names(dds_CA1_list), 'All')) {
  print(tissue)
  #We extract the respective Deseq2 object
  dds_tissue <- dds_CA1_list[[tissue]]
  
  #We use the DESE2-normalized data, which we extract by getting the counts from the Deseq2 object
  count_table <- counts(dds_tissue, normalized=T)
  
  count_table <- as.data.frame(count_table)
  
  
  #To make sure that we are not correlating flat noise with flat noise, we will focus on genes that
  #we consider as 'expressed', which means we'll discard all genes that Deseq2 has flagged in their independent filtering by assigning a na in the padj column
  #We'll use the results list from 3 to 28 months to make sure we capture genes only expressed at old age
  ExpresedGenes <- results_list_CA1$Both[[tissue]]$`28_vs_3`$resOrdered %>% filter(!is.na(padj)) %>% pull(gene_symbol)
  #We definie the numberic values we want to correlate against - in this case, the age of the mice
  Score_to_compare <- as.numeric(as.character(as.data.frame(colData(dds_tissue))[,'age']))
  
  #Perform gene-vs-age score correlataion
  #To do that we set up a dataframe that can holde the results
  #We want to record both the spearman and pearson correlation corefficient, as well as the results from cor.test to see that we have statistical
  #confidence in our correlation coefficents
  correlation_frame <- data_frame(gene='A', cor_spear=1, p_spear=1, cor_pears=1, p_pears=1)
  
  #We set up a counter so that we can have a progressbar for our analysis
  i <- 1
  pb = txtProgressBar(min = 1, max = length(ExpresedGenes), initial = 1, style = 3)
  #We will now iterate over very expresed gene and perform the following steps
  for (gene in ExpresedGenes) {
    setTxtProgressBar(pb,i)
    #We set up a temporary dataframe that holds age and normalized counts 
    cor_frame <- data.frame(expr=t(count_table[gene,]), concentration= Score_to_compare)
    colnames(cor_frame) <- c('expr', 'concentration')
    #We wil only retain samples where the count is not 0
    cor_frame <- cor_frame %>% filter(expr > 0)
    #After filtering for non-expressing samples, we will only continue the analysis if we have at least 20 samples letf - otherwise we end up
    #with very spurious correlations
    if (nrow(cor_frame) >= 20) {
      #Now we run the correlation tests for spearman and pearson
      cor_vector_spearman <- cor.test(cor_frame$expr, cor_frame$concentration, method = 'spearman')
      cor_vector_pearson <- cor.test(cor_frame$expr, cor_frame$concentration, method = 'pearson')
      #and we now appending to the correlation_frame a new row that contains all the coeffiients and pvalues
      correlation_frame[i, 'gene'] <- gene
      correlation_frame[i, 'cellsExpr'] <- nrow(cor_frame)
      correlation_frame[i, 'cor_spear'] <- cor_vector_spearman$estimate
      correlation_frame[i, 'p_spear'] <- cor_vector_spearman$p.value
      correlation_frame[i, 'cor_pears'] <- cor_vector_pearson$estimate
      correlation_frame[i, 'p_pears'] <- cor_vector_pearson$p.value

      #if the gene had fewer than 20 samples expressing it, we'll add a dummy row
    } else {
      correlation_frame[i, 'gene'] <- gene
      correlation_frame[i, 'cellsExpr'] <- nrow(cor_frame)
      correlation_frame[i, 'cor_spear'] <- NA
      correlation_frame[i, 'p_spear'] <- NA
      correlation_frame[i, 'cor_pears'] <- NA
      correlation_frame[i, 'p_pears'] <- NA
      
    }
    #We add 1 to the conter for the progress bar
    i <- i + 1
  }
  #Cleanup results and run the multiple testing correction
  correlation_frame <- na.omit(correlation_frame)
  correlation_frame$padj_spear <- p.adjust(correlation_frame$p_spear, method = 'BH')
  correlation_frame$padj_pear <- p.adjust(correlation_frame$p_pears, method = 'BH')
  #Relabel
  tissue_correlation_fullTable <- correlation_frame
  rm(correlation_frame)
  #Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibhitng an absolute correlation value of 0.5)
  print(paste('TopCorGenes: ', nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
  nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05))
  #Store data the results in our correlation results list
  BulkSeq_Agecorrelation_results_list[[tissue]][['resall']] <- tissue_correlation_fullTable
  #Filter and store data for high correlates in our correlation results list
  BulkSeq_Agecorrelation_results_list[[tissue]][['topCor']]  <- (tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))


  
}
#We will save the results list for later usage
save(BulkSeq_Agecorrelation_results_list, file='R_objects/dds_BulkSeq_Aging_CorrelationResults.bin')


#Now that we've run the correlation analysis we want to plot how many genes correlated significanlty with age in each tissue/region

#We will set a range of pvalue cutoffs, to see if the results remain stable 
padj_range <-c(0.001, 0.01, 0.05, 0.1)
#We have to set a correlation coefficient cutoff (absolute cutoff of 0.5)
cor_cutoff <- 0.5

#We will have to establish a 'master data frame' with some dummy data. We will iteratively add data to this data frame
plot_df_crossTissue <- data.frame(padj_cutoff=0.05, Type='A', Comparison='A', Numb_of_genes=1, tissue='A')

results_list_CA1_bothSex <- BulkSeq_Agecorrelation_results_list
#We will iterate over each of the tissues, get the correlation results to build the bargraph plots
for (tissue in names(results_list_CA1_bothSex)) {
  
  print(paste('Processing tissue:', tissue, sep = ' '))
  #Extract the relevant resultstable
  results_list_tissue <- results_list_CA1_bothSex[[tissue]]
  #We'll have to set up a list onto which we'll attach the extracted number of DEGs to build the bargraph plots
  barplot_list <- list()
  #We'll repeat this step for each p value cutoff
  for (padj_cutoff in c(0.001,0.01,0.05,0.1)) {
    barplot_frame <- data.frame()
    comparison <- 'correlation'
    #We extract the significantly postiively and negatively correlated genes seperately and store their total number in the bargraph dataframe
    results_list_tissue_temp <- as.data.frame(results_list_tissue$resall)
    down_regulated <- nrow(results_list_tissue_temp %>% filter(padj_spear < padj_cutoff) %>% filter(cor_spear < -1*cor_cutoff))
    up_regulated <- nrow(results_list_tissue_temp %>% filter(padj_spear < padj_cutoff) %>% filter(cor_spear > 1*cor_cutoff))
    barplot_frame["down",comparison] <- data.frame(down_regulated)
    barplot_frame["up",comparison] <- data.frame(up_regulated)
    #After having the respective information assembled, we'sll store the outcome
    barplot_list[[as.character(padj_cutoff)]] <- barplot_frame
    
  }
  #Having assemled the relevant information across a range of pvalue cutoffs, we can now aggregate the results
  barplot_list <- bind_rows(barplot_list,.id = "id")
  #We'll use simple labels of to define the direction. Since we strictly followed the design of getting first down/up-regulated genes, we'll just apply the rep function
  barplot_list$Type <- rep(c('down','up'), length(padj_range))
  #Time to 'melt' the dataframe into the long format
  plot_df <- melt(barplot_list, id.vars = c('id','Type'))
  #Relabel columns
  colnames(plot_df) <- c("padj_cutoff", 'Type',"Comparison", "Numb_of_genes")
  
  plot_df$tissue <- tissue
  #Having gathered and prepared all the relevant information from this tissue/region, we're now appending that to the dataframe created previously
  plot_df_crossTissue <- rbind(plot_df_crossTissue, plot_df)
  
  
}

plot_df_crossTissue <- plot_df_crossTissue[-1,] #Remove the first row we needed to set up the dataframe

#we'll also have to set the order of the factors in the 'Type' column that contains informatino about  up-/down-regualtion
plot_df_crossTissue$Type <- factor(plot_df_crossTissue$Type, levels = c("up","down"), ordered = T)

#Now order the tissues according to their total number of age-correlated genes
plot_df_crossTissue$tissue <- factor(plot_df_crossTissue$tissue, levels = rev(as.character((plot_df_crossTissue  %>% dplyr::filter(padj_cutoff == 0.05) %>% dplyr::group_by(tissue, padj_cutoff)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
                                     ordered = T)
#For the resulting plot, we'll focus on the genes passing the padj cutoff of 0.05 
plot_df_crossTissue <- plot_df_crossTissue %>% filter(padj_cutoff == 0.05)
#####And plot Figure 1H
ggplot(plot_df_crossTissue, aes(x=tissue, y=Numb_of_genes, fill=Type), color=NA) + geom_bar(stat = "identity") 

