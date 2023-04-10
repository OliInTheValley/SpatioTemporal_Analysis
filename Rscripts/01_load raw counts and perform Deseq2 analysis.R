rm(list=ls())


library(readr)
library(DESeq2)

###Start
#Load dataset from raw counttable
library(readr)
BulkSeq_Aging_Counttable <- read_delim("input_data/BulkSeq_Aging/BulkSeq_Aging_Counttable.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)
BulkSeq_Aging_Counttable <- as.data.frame(BulkSeq_Aging_Counttable)
row.names(BulkSeq_Aging_Counttable) <- BulkSeq_Aging_Counttable[,1] #turn gene names into rownames of the table
BulkSeq_Aging_Counttable[1] <- NULL

#Load meta data and create an experimentDesign frame that fits Deseq2's design
experimentDesign  <- data.frame(row.names = colnames(BulkSeq_Aging_Counttable)) #Create dataframe with sampleIDs as rownames

BulkSeq_Aging_metadata <- read_delim("input_data/BulkSeq_Aging/BulkSeq_Aging_metadata.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
BulkSeq_Aging_metadata <- as.data.frame(BulkSeq_Aging_metadata)

#Transfer information from the meta data table into a Deseq2-compatible experimentDesign structure
for (meta_column in colnames(BulkSeq_Aging_metadata)[-which(colnames(BulkSeq_Aging_metadata)=='sampleID')]) {
  experimentDesign[,meta_column] <- plyr::mapvalues(x = row.names(experimentDesign), from = BulkSeq_Aging_metadata[,'sampleID'], to = BulkSeq_Aging_metadata[,meta_column])
}

(experimentDesign <- as.data.frame(unclass(experimentDesign), row.names = row.names(experimentDesign)))

#For simplicity, we will create a column called 'tissue' which contains acronyms for the verbose tissue names currently provided in the metadata

long_labels <- c("Corpus callosum",
                 "Cerebellum",
                 "Motor cortex",
                 "Caudate putamen",
                 "Entorhinal cortex",
                 "Hippocampus (posterior)",
                 "Hippocampus (anterior)", 
                 "Hypothalamus", 
                 "Medulla",
                 "Olfactory bulb", 
                 "Choroid Plexus",
                 "Pons",
                 "Subventricular zone",
                 "Thalamus",
                 "Visual Cortex")

acronym_labels <- c("cc",
                    "cer",
                    "cor",
                    "cp",
                    "ent",
                    "hi2",
                    "hi", 
                    "hy", 
                    "med",
                    "olf", 
                    "plx",
                    "pon",
                    "svz",
                    "th",
                    "vis")
#Relabel 
experimentDesign$tissue <- plyr::mapvalues(x = experimentDesign$tissueLong, from = long_labels, to = acronym_labels)


#Create Deseq2 object
#This will be the main processed data object to work from. We use the acronym 'CA' for "cerebrum aevum" (latin for aging brain)
#to label the samples from the aging bulkseq dataset, we use numer '1', hence CA1
dds_CA1 <- DESeqDataSetFromMatrix(countData=BulkSeq_Aging_Counttable, colData=experimentDesign, design = ~age + tissue ) 
#we first run the analysis using a simple model without an interaction term
design(dds_CA1) <-  ~age + tissue
#Run Deseq2
dds_CA1 <- DESeq(dds_CA1, fitType = 'local')
#This will create a list onto which we will add reiteratively objects for each/region individually
dds_CA1_list <- list(All=dds_CA1
)


#Now we create a results list onto which we will reiterativley add results tables from the Deseq2 analysis
results_list_CA1 <- list()
results_list_CA1$Both <- list()

#Define cutoff for differntial expression test 
padj_cutoff <- 0.05
#Set the 'sex' that hsould be analysed in this run (choose 'Both', 'Male', 'Female')
sex <- 'Both'

#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective region and then run the regular Deseq2 workflow
for (tissue in c( unique(as.character(dds_CA1_list$All$tissue)))) {
  print(tissue)
  #subset for current tissue and re-run Deseq2
  if (tissue != 'All') {
    dds_CA1_list[[tissue]] <- dds_CA1[, dds_CA1$tissue == tissue]
    dds_CA1_list[[tissue]]$tissue <- droplevels(dds_CA1_list[[tissue]]$tissue)
    design(dds_CA1_list[[tissue]]) <- ~age + sex
    #Run Deseq2 and store in the Deseq2 object list
    dds_CA1_list[[tissue]] <- DESeq(dds_CA1_list[[tissue]], fitType = 'local')
  }
  
  #Now take the Deseq2 object and run the extraction of the pairwise comparisons
  dds_tissue_temp <- dds_CA1_list[[tissue]]
  results_list_tissue <- list()
  
  #Create a list of all pairwise comparisons
  comparison_list <- t(combn(unique(colData(dds_tissue_temp)$age),2))
  
  #Iterate over all pairwise comparisons and extract the resultstable, store it in the results list. then flip the conditions (so that results for 3 vs 21 and 21 vs 3 months is stored)
  for (row in 1:nrow(comparison_list)) {
    print(row)
    #Get the conditions/timepoints to test
    cond1 <- as.character(comparison_list[row,1])
    cond2 <- as.character(comparison_list[row,2])

    folder <- paste(cond1, cond2, sep = "_vs_")
    print(folder)
    
    aspect_to_test <- 'age'
    
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
    
    aspect_to_test <- 'age'
    
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
  results_list_CA1[[sex]][[tissue]]  <- results_list_tissue
  #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
  save(dds_CA1_list, results_list_CA1,  file='R_objects/dds_BulkSeq_Aging.bin')
  
}



#Now that we've run Deseq2 across all age groups witin each individual tissue, we'll explore the number of significant DEGs per pair-wise comparison, per tissue
#We will limit our diagnostic plots to genes that are flagged as significant in at least 2 pairwise comparisons from 3 vs any older age (e.g. 3 vs 12 and 3 vs 21)

#We will set a range of pvalue cutoffs, to see if the results remain stable 
padj_range <-c(0.001, 0.01, 0.05, 0.1)

#We will have to establish a 'master data frame' with some dummy data. We will iteratively add data to this data frame
plot_df_crossTissue <- data.frame(padj_cutoff=0.05, Type='A', Comparison='A', Numb_of_genes=1, tissue='A')

#We will focus here only the results whwere we combined both sexes
results_list_CA1_bothSex <- results_list_CA1$Both
comparisons_to_plot <- c('12_vs_3', '15_vs_3',
                         '18_vs_3', '21_vs_3', '26_vs_3',
                         '28_vs_3')
for (tissue in names(results_list_CA1_bothSex)) {
  
  print(paste('Processing tissue:', tissue, sep = ' '))
  
  #First we'll get the current tissue's results tables from our results list
  results_list_tissue <- results_list_CA1_bothSex[[tissue]]
  
  #We'll have to set up a list onto which we'll attach the extracted number of DEGs to build the bargraph plots
  barplot_list <- list()
  #We'll repeat this step for each p value cutoff
  for (padj_cutoff in padj_range) {
    
    barplot_frame <- data.frame()
    
    for (comparison in comparisons_to_plot) {
      #First, we'll have to make sure that we only focus on genes that have found at least one other comparison
      #To that end, we'll get DEGs from all other comparions except the one we're currently investigating, which is exactly what our everything_but() function does
      other_comparisons <- everything_but(c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3','28_vs_3'), comparison)
      #Let's extract the DEGs (as defined by the curren pvalue cutoff) from the comparison we want to examine
      #We will do this first for all the down-regulated genes, and later for all the up-regualted genes
      results_list_tissue_temp <- as.data.frame(results_list_tissue[[comparison]]$resall)
      #lets take the down-regulated genes from the comparison we're intersted in as reference, or 'target'
      down_regulated_genes_target <- (results_list_tissue_temp %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange<0) %>% pull(gene_symbol))
      #and now we go over every other comparison but the one we're looking at and get every possible DEG in the same direction from every available pair-wise comparison
      down_regulated_genes_background <- c()
      for (background_comparsion in other_comparisons) {
        results_list_tissue_background <- as.data.frame(results_list_tissue[[background_comparsion]]$resall)
        down_regulated_genes_background <- c(down_regulated_genes_background, (results_list_tissue_background %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange<0) %>% pull(gene_symbol)))
      }
      #Our 'background' is in this case all the DEGs from pairwise comparisons other the one we're currently looking at
      #We have assembled those in the previous loop. We only need every DEG named once, so we'll turn our vector into one with unique elements
      down_regulated_genes_background <- unique(down_regulated_genes_background)
      #Now comes the crux. DEG of our current comaprions is only a 'true' DEG by our definition if it has also found in at least one of the other 'background' comparisons
      #Only genes present in our current vector and the background vector will be retained. That makes our down_regulatd genes
      #Here, we're only intersted in how many DEGs are retained, so we'll just extract the length of that intersect
      down_regulated <- length(intersect(down_regulated_genes_target, down_regulated_genes_background))
      
      #Okay that worked well, now we'll repeat that for the genes that have a positive log2Foldchange
      up_regulated_genes_target <- (results_list_tissue_temp %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange>0) %>% pull(gene_symbol))
      up_regulated_genes_background <- c()
      for (background_comparsion in other_comparisons) {
        results_list_tissue_background <- as.data.frame(results_list_tissue[[background_comparsion]]$resall)
        up_regulated_genes_background <- c(up_regulated_genes_background, (results_list_tissue_background %>% filter(padj < padj_cutoff) %>% filter(log2FoldChange>0) %>% pull(gene_symbol)))
      }
      up_regulated_genes_background <- unique(up_regulated_genes_background)
      #Same as before,we're only intersted in how many DEGs are retained, so we'll just extract the length of that intersect
      up_regulated <- length(intersect(up_regulated_genes_target, up_regulated_genes_background))
      
      #Now we'll report the number of up/down-regulated DEGs into a dataframe, using the comparison as a new column
      barplot_frame["down",comparison] <- data.frame(down_regulated)
      barplot_frame["up",comparison] <- data.frame(up_regulated)
      #The dataframe we generated is only relevant for the current value cutoff - so we'll add that to the list we generated before starting this loop
      barplot_list[[as.character(padj_cutoff)]] <- barplot_frame
    }
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

#we'll have to define the order of the comparisons, to make sure the resulting plot arranges the x axis properly
comparison_order <- c('12_vs_3', '15_vs_3',
                         '18_vs_3', '21_vs_3', '26_vs_3',
                         '28_vs_3')
plot_df_crossTissue$Comparison <- factor(plot_df_crossTissue$Comparison, levels=comparison_order,
                                         ordered=T)
#we'll also have to set the order of the factors in the 'Type' column that contains informatino about  up-/down-regualtion
plot_df_crossTissue$Type <- factor(plot_df_crossTissue$Type, levels = c("up","down"), ordered = T)



##Now we'll set up the data to plot a line graph of up/down-regualted DEGs 
#For this we'll focus only on DEGs passing the cutoff of 0.05
plot_df_crossTissue_line <- plot_df_crossTissue %>% filter(padj_cutoff == 0.05) 
#We'll give the number of genes that account for the down-regulated genes a negative sign to disentangle them later in the line plot
plot_df_crossTissue_line[which(plot_df_crossTissue_line$Type=='down'),'Numb_of_genes'] <- plot_df_crossTissue_line[which(plot_df_crossTissue_line$Type=='down'),'Numb_of_genes']*-1
#We now need some grouping info. the easiest is to make a new column containing a combined strin of tissue and 'type' (ie regulation)
plot_df_crossTissue_line$type_combined <- paste(plot_df_crossTissue_line$tissue, plot_df_crossTissue_line$Type, sep = '_')

#We only have 6 'age point's (7 age groups minus 1, because we reference everything to 3 months)
#This is too few for ggplots smoothing function if a categorial variable (in this case the comparison eg. "12 vs 3 months") is used. So we'll create dummy values from 1 to 6 for this
lat.order.levels <- comparisons_to_plot
#we provide a dataframe where we indicate which of the continious dummy values correspond to which comparison
site.data <- data.frame(
  Comparison = lat.order.levels,
  site_num = 1:length(lat.order.levels)
)
#We'll apply the same factor order to this dataframe
site.data$Comparison <- factor(site.data$Comparison, levels=comparison_order,
                               ordered=T)

#Now we'll join these
df <- dplyr::left_join(plot_df_crossTissue_line, site.data, by= "Comparison")
#####And plot Figure 1F
(myplot <- ggplot(df, aes(x=site_num, y=Numb_of_genes, colour=tissue, group=type_combined)) +
    geom_smooth(se = F) +
    scale_x_continuous(
      breaks=1:length(lat.order.levels),
      labels=lat.order.levels, expand = c(0,0)
    ) +  scale_y_continuous(expand = c(0,0), limits = c(-1000,1000),name='# of diff. genes') 
)

#If we want to plot this information as heatmap we do the following
#First we split above's dataframe into down regulation information only
plot_df_crossTissue_line_down <- plot_df_crossTissue_line %>% filter(Type == 'down')
#Now we reformat that using tidyr's spread function to go from long to wide format
plot_df_crossTissue_line_down_wide <- spread(plot_df_crossTissue_line_down, Comparison, Numb_of_genes)
#We will need only the 3rd column (tissue) and the information per comaprions (5th to 10th column)
plot_df_crossTissue_line_down_wide <- plot_df_crossTissue_line_down_wide[c(3,5:10)]
#Assigning the tissue information as rownames for easier plotting
row.names(plot_df_crossTissue_line_down_wide) <- plot_df_crossTissue_line_down_wide$tissue

#Now we repate that for the up-regualted genes
plot_df_crossTissue_line_up <- plot_df_crossTissue_line %>% filter(Type == 'up')
plot_df_crossTissue_line_up_wide <- spread(plot_df_crossTissue_line_up, Comparison, Numb_of_genes)
plot_df_crossTissue_line_up_wide <- plot_df_crossTissue_line_up_wide[c(3,5:10)]
row.names(plot_df_crossTissue_line_up_wide) <- plot_df_crossTissue_line_up_wide$tissue

#Now we'll cbind the two dataframes
plot_df_crossTissue_line_comb_wide <- cbind(plot_df_crossTissue_line_down_wide[c(comparisons_to_plot)], plot_df_crossTissue_line_up_wide[c(comparisons_to_plot)])


#####And plot Figure 1G
pheatmap(plot_df_crossTissue_line_comb_wide, 
         show_colnames = T, clustering_method = 'ward.D2')




#We will also create a bargraph that summarises all DEGs per regions (passing the padj cutoff in at least two comparisons with 3m)
#We set a padj cutoff of 0.05 to mark a gene as differentially expressed.
padj_cutoff <- 0.05

plot_df_crossTissue <- data.frame(direction='A', Comparison='A', Numb_of_genes=1, tissue='A', sex='A')
#We will focus here on the data that covers both sexes
results_list_CA1_bothSex <- results_list_CA1$Both
#We'll iterate over each tissue, collapse the results lists from all the age-related comparisons and count in how many of these 
#a given gene passed the significance threshold

for (tissue in names(results_list_CA1_bothSex)) {
  
  print(paste('Processing tissue:', tissue, sep = ' '))
  #We subset for the results tables from the current tissue
  results_list_tissue <- results_list_CA1_bothSex[[tissue]]
  #We create a dataframe that will store the information 
  barplot_frame <- data.frame()
  
  #We will need to know in which direction a given gene is moving (up/downregualted)
  #To this end, we'll take the 21 vs 3 months-comaprison as refernce and use the sign of the log2FC as indicator of the regulation dirextion
  regulation_df <- results_list_tissue[["21_vs_3"]]$resOrdered
  #We use lappy to collapse the information from all 6 age-related comparisons
  DEGs_of_tissue <- (unlist((lapply(results_list_tissue[c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3', '28_vs_3' )], function (x) (row.names(subset(x[["resall"]], padj < padj_cutoff)))))))
  #now we count often a given gene 'occured'  in the DEGS_of_tissue vector and supset for genes that have been found at least twice (i.e. passing two times the significance threshould)
  DEGs_of_tissue <- names(which(table(DEGs_of_tissue) >= 2))
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

#Now we plot this and generate Figure S4A
ggplot(plot_df_crossTissue , aes(x=tissue, y=Numb_of_genes, fill=direction), color=NA) + geom_bar(stat = "identity")

