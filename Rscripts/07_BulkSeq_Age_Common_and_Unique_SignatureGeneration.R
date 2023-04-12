
###Here, we want to inspect which DEGs are shared between tissues/regions and which are unique
#We then build signatures for common and region-specific genes
#We load the DEseq2 object containing the bulk seq results
load('R_objects/dds_BulkSeq_Aging.bin')




#First we assemble in which tissue which gene is passing the definition of an age-related DEG (passing sinificance threshold in at least two comparisons to 3m)
#our method of getting an overview is the upset plot. To create one, we will need a matrix of tissues (columns) by genes (rows) where the cell indicate whether the given gene
#is a DEG in that tissue (DEG is defined to pass two different age-comarisons in a given tissue)
#We set the significance threshold
padj_cutoff <- 0.05
#And the comparisons we want to use to define the DEGs
comparisons_to_analyze <-  c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3', '28_vs_3')
#We set the number comparisons a gene needs to pass statistical testing 
min_Difftimepoints <- 2
#We focus only on the results coming from analysing both sexes together
results_list_CA1_bothSex <- results_list_CA1$Both  
#We will need to establish a list that we'll add the extracted information; here we'll add a table for each tissue where we summarise for each gene informatino gathered across all 6 age-realted comparisons
CA1_AgingDEGS_joined_list <- list()
#We set up a dummy dataframe onto which we will append results
comb_data_frame <- data.frame('DG'=1, row.names = 'SomeGene' )
#We exclude the hippocampus (posterior) from this analysis so that we won't be biased in favor or against hippocampus (as we have two tissue punches from that region)
#For this purpose we exclude it using the everything_but function
for (tissue in everything_but(names(results_list_CA1_bothSex), c('hi2'))) {
  print(tissue)
  #In order for the later plotting to work, we will need tables with the exact same genes (i.e. all quanitfied ones) with some information which one was an age-related DEG
  #We will iterate over all comparisons (age vs 3 months) and assemble the resultstables and later aggregate the information
  #We will need a dummy dataframe to append to. Therefore we extract first comparison 12 vs 3 and build basic dataframe and append
  #Only take genes with a padj beneath the cutoff to ensure no Deseq2 "NA" data is introduced 
  comparison_df <- head(subset(as.data.frame(results_list_CA1_bothSex[[tissue]][[comparisons_to_analyze[1]]]$resall), padj <= 1), n=1)
  #We need a column in which we write what comparison this gene has been coming from. For this first line, we use a dummy value "init"
  comparison_df$tested_comparison <- 'init'
  for (comparison in comparisons_to_analyze) {
    #Well get the resultstable of that particular comparison and append
    comparison_df_temp <- as.data.frame(as.data.frame(results_list_CA1_bothSex[[tissue]][[comparison]]$resall))
    comparison_df_temp$tested_comparison <- comparison
    comparison_df <- rbind(comparison_df, comparison_df_temp)
  }
  #First we remove the initial row
  comparison_df <- comparison_df %>% filter(tested_comparison != 'init')
  #Then we use dplyr to aggregate the information regarding in how many comparisons a gene has been crossing the significance threshould ('occured')
  #We'll also gather ifnormation about the median log2FC across all comparisons, as well as the baseMean. As recoreded padj, we will the lowest one of all significant comarsons
  temp_df <- comparison_df %>% 
    dplyr::group_by(gene_symbol) %>% 
    dplyr::summarise(occurance = sum(padj <= padj_cutoff), log2FoldChange=median(log2FoldChange), padj=min(padj) ,baseMean=median(baseMean))
  #Assign if a gene is an age-related DEG if it passed the minimally requied number of comparisons the significance cutoff
  temp_df <-  temp_df %>% mutate(AgeRelated_DEG=ifelse(is.na(padj), 0, ifelse(occurance >=  min_Difftimepoints, ifelse(padj < padj_cutoff, 1, 0), 0) ))   
  #Likewise, we'll use the log2FC information to define if the gene is up/down-regualted or not affected (nm)
  temp_df$direction <- unlist(lapply(temp_df$log2FoldChange, function(x) ifelse(is.na(x), 'nm', ifelse(x > 0, 'up', 'down'))))
  temp_df <- as.data.frame(temp_df)
  #Well set the genes as row.names to fascilitate the later reformatting
  row.names(temp_df) <- temp_df$gene_symbol
  #We select and carry forward the columns that are critical to us
  temp_df <- temp_df[c('AgeRelated_DEG', 'direction', 'occurance', 'gene_symbol', 'log2FoldChange', 'padj', 'baseMean')]
  #The result should be a dataframe that has always the same length (21,076 genes), ie all quantified genes. This is essential to be the case across all tissues
  CA1_AgingDEGS_joined_list[[tissue]] <- temp_df
  
}
#we use cbind.data.frame to collapse all dataframes from each tissue together
plot_df <- cbind.data.frame(CA1_AgingDEGS_joined_list)

#Now we'll extract all the columns that contain the ifnromation about wether the gene is an age-related DEG
#The reuslt will be a matrix that contains the tissues as columns and the rows corresponding to a gene, with values of 0/1 indicating if a gene is a DEG
plot_df_cond <- plot_df[,grepl('AgeRelated_DEG',colnames(plot_df))]
#We'll have to reformat a bit
colnames(plot_df_cond) <- gsub('\\.AgeRelated_DEG', '', colnames(plot_df_cond))
#If we have done everything correctly, then all that's left afer the reformatting should be the 15 regions 
tissues_to_plot <- colnames(plot_df_cond)

######The result should look like this:
#               cc cer cor cp ent hi hy med olf plx pon svz th
# 1700049E17Rik2  0   0   0  0   0  0  0   0   0   0   0   0  0
# 9330162012Rik   0   0   0  0   0  0  0   0   0   0   0   0  0
# a               0   0   0  0   0  0  0   0   0   0   0   0  0
# A1bg            0   0   0  0   0  0  0   0   0   0   0   0  0
# A1cf            0   0   0  0   0  0  0   0   0   0   0   0  0
# A26c2           0   0   0  0   0  0  0   0   0   0   0   0  0
# A2m             0   0   0  0   0  0  0   1   0   1   0   0  0
# A2ml1           0   0   0  0   0  0  0   0   0   0   0   0  0
# A3galt2         0   0   0  0   0  0  0   0   0   0   0   0  0
# A4galt          0   0   0  0   0  0  0   0   0   1   0   0  0

#Now we have a convenient format that allows direct comparisosn of DEG-detection across 15 regions
#First we'll subset this dataframe to genes that are at in at least 1 tissue an age-related DEG
plot_df_cond_min <- plot_df_cond[which(rowSums(plot_df_cond[tissues_to_plot])>0),]
plot_df_cond_min <- as.data.frame(plot_df_cond_min)
#This is a very useful representation of all the DEG results so we'll label it properly and store it as an R object
CA1_AgingDEGS_genetable_min <- plot_df_cond_min

#The joind_list that contains for each tissue a single table that summarises the results from across all age-realted comparisons will be useful in the future. We'll save it together  
#With the DEG_genetable_min

save(CA1_AgingDEGS_genetable_min, CA1_AgingDEGS_joined_list, file='R_objects/BulkSeq_AgingDEGS_genetable.bin')

#There should be only rows left where at least one tissue exhibits a 1 (ie age-related DEG)
#We create a dataframe with boolean values
plot_df_cond_min_boolean <- CA1_AgingDEGS_genetable_min[tissues_to_plot] == 1

#The output should look like this
#                cc   cer   cor    cp   ent    hi    hy   med   olf   plx   pon   svz    th
# A2m           FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE
# A4galt        FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
# AA414992      FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
# AA467197      FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE FALSE FALSE
# Aacs          FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Aagab         FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# Aars          FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

#We will now inspect in how many tissues/regions a given gene was detected as an age-related DEG
#We use rowSums and sort the resutls
gene_count <- sort(rowSums(CA1_AgingDEGS_genetable_min), decreasing = T)
#We turn this into a dataframe and plot this as a bargarph
gene_count_df <- data.frame(count=gene_count, row.names = names(gene_count))
#Now we plot this and generate Figure 2A
ggplot(data=gene_count_df, aes(x=count)) + geom_bar() + scale_x_continuous(limits = c(0,15), expand = c(0,0), breaks = seq(1,15,1))

#We can use the boolean table as input for the upset plot function
#We will limit the output to sets that have at least 25 genes
#This creates Fig 6A - we will need this for the later analysis of region-specific signatures
ComplexUpset::upset(plot_df_cond_min,tissues_to_plot, width_ratio=0.1,min_size=25,
                    base_annotations=list(
                      'Intersection size'=intersection_size(
                        counts=T 
                      )
                    )
)



### Now that we have assembled the required information across tissues, we can build our signatures
#First, let's generate a signature that captures genes that are found across 10 regions min - this will form our Common Aging Signature (CAS). We will use VISION once again for this
#We define and store our CAS genes, as we'll need them multiple times in the future
CA1_CASgenes <- names(which(rowSums(CA1_AgingDEGS_genetable_min) >= 10))
save(CA1_CASgenes, file='R_objects/CA1_CASgenes.bin')
#First, we'll set up the signature vector
mySignatures_CAS_AgingDEGs_perRegion <- c()
#This vector will cotain both up/down-regulated genes; to get the direction of the CAS genes, we'll use the corpus callosum (as these genes are all regulated in the same direction 
#across tissues/regions), any other region could also be used; we'll take the 26 vs 3 months comparison as represntative

#VISION requires a named vector where the values are -1 or 1 (for up/down regulation); so we'll create a vector with +1 of the same length as the marker genes
signature_vector <-  results_list_CA1$Both$cc$`26_vs_3`$resOrdered[CA1_CASgenes,'log2FoldChange']
#We'll transform the vector into a simple -1/+1 vector
signature_vector[which(signature_vector<0)] <- -1
signature_vector[which(signature_vector>0)] <- 1
#Now we'll name that vector. Vision requirs all genes to be uppercase
names(signature_vector) <- toupper(CA1_CASgenes)
#We'll remove any potential dupplicates that could cause trouble down the line
signature_vector <- signature_vector[!(duplicated(names(signature_vector)) | duplicated(names(signature_vector),fromLast=TRUE))]

sig_tissue_markergenes <- createGeneSignature(name = 'CommonAgingSignature', sigData = signature_vector)
#We will append the signature to the signature vector
mySignatures_CAS_AgingDEGs_perRegion <- c(mySignatures_CAS_AgingDEGs_perRegion, sig_tissue_markergenes)


##Let's visualize some of the CAS genes as a heatmap
genes_to_plot <- CA1_CASgenes
#We'd like to show that these genes incarese/decrease with age over time. So we'll extract for each tissue and each of the CAS genes the log2FoldChange 'sequence' from young to old
#Ie we  get the log2FC from 12 vs 3, 15 vs 3, 18 vs 3 etc
#We'll use a dataframe to store all the relevant information - imporant is that it contains the comparisons as rownames
plot_masterFrame <- data.frame(init=rep(0, length(genes_to_plot)), row.names = genes_to_plot)
for (tissue in names(results_list_CA1[['Both']])) {
  #We'll subset the results list for the respective tissue and we'll focus on the data coming from comparisons where male/female are pooled
  results_list_temp <- results_list_CA1[['Both']][[tissue]]
  results_list_temp <- results_list_temp[c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3', '28_vs_3')]
  #We'll lapply to assemble for all of the genes the log2FCs as a single dataframe, where columns are the genes, rows are the comparisons
  input_table <- do.call(rbind, (lapply(results_list_temp, function (x) x[["resall"]][genes_to_plot,'log2FoldChange'])))
  #We'll relable the data frame
  colnames(input_table) <- c(genes_to_plot)
  row.names(input_table) <- gsub('_vs_3', '', row.names(input_table))
  #Finally, we tranform the dataframe so that 'age' (vs 3months) demarcates the columns
  input_table <- as.data.frame(t(input_table))
  #We will add the (current) tissue name to the 'age' label so we can use that later to arrange the heatmap
  colnames(input_table) <- paste(colnames(input_table), tissue, sep = '_')
  #If we did everything correctly this is what the output could/should look like (for the visual cortex)
  #             12_vis       15_vis       18_vis       21_vis      26_vis       28_vis
  # Abca8a    3.722765e-01  0.410002166  0.598466999  0.720168232  0.44963548  1.026053152
  # Ahsa1    -4.596715e-02 -0.034742620 -0.044496188 -0.062060912 -0.03390148 -0.156564525
  # Aif1      3.116152e-01  0.327878417  0.306134376  0.412190306  0.53909068  0.136117476
  # B2m       6.775408e-02  0.216166081  0.098467926  0.288375223  0.26961020  0.241254402
  # Bcl2a1b   3.986364e-01  0.596961686  0.297591443  0.169210353  0.80591939  0.64175758
  
  plot_masterFrame <- cbind(plot_masterFrame, input_table)
}
#We remove the inital 'dummy' column
plot_masterFrame$init <- NULL
#Time to plot the heatmap
#First we set a color palette
colfunc <- colorRampPalette(c("#19F5FF","black", "#ffff36"))
#We have to set here the min/max value
min_value <- -1.5
max_value <- 3

color_frame <- c(min_value, max_value)
#Now we're building our palaette and arrange it so that it represents the overall picture properly
paletteLength <- 100
myColor <- colorRampPalette(c("#19F5FF","black", "#ffff36"))(paletteLength)

myBreaks <- c(seq(min(color_frame), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(color_frame)/paletteLength, max(color_frame), length.out=floor(paletteLength/2)))
#Let's plot the heatmap, which results in Figure 2B
pheatmap(plot_masterFrame, 
         show_colnames = T, 
         cluster_cols = F,  
         color= myColor, breaks=myBreaks, clustering_method = "ward.D2", border_color = NA)
#We can also subset the dataframe to fewer genes, like the ones that made it past the threshold in 12 instead of 10 gene
pheatmap(plot_masterFrame[names(which(rowSums(CA1_AgingDEGS_genetable_min) >= 11)),] , 
         show_colnames = T, 
         cluster_cols = F,  
         color= myColor, breaks=myBreaks, clustering_method = "ward.D2", border_color = NA)


###So much for the common genes. We can use the sama material that got us the CAS to define region-specific signatures
#Essentially, we'll extract the same information - the number of unique DEGs for each region - what marks the foundation of the UpSet plot
#So we'll select the unique DEGs and create signatures for each tissue

#First, we'll define the 'unique' gene set. We ask which genes have a rowsum in the AgingDEGs_genetable with only a 1
unique_DEGs <-  names(which(rowSums(CA1_AgingDEGS_genetable_min) <= 1))


#Now we'll rotate over each tissue, and build a signature, given that the number of DEGs in the signature is passing a threshold of 20
#We will do this for signatures based on up- and dowwn-regulated genes, only up-regulated, or only down-regulated genes. Always with the cutoff of at least 20
#We will use the joined-list for this puprose
for (tissue in names(CA1_AgingDEGS_joined_list)) {
  #As always with VISION we need to supply which direcitno the gene is going (up/down) as demarcated by a +1/-1, respectively
  #We subset for the respective tissu
  results_df <- as.data.frame(CA1_AgingDEGS_joined_list[[tissue]])
  #Only consider genes that have been defined as 'AgeRelated_DEG'
  results_df <- results_df %>% dplyr::filter(AgeRelated_DEG == 1)
  #Filter for unique genes
  results_df <- results_df %>% filter(gene_symbol %in% unique_DEGs)
  #Now we'll assemble the signature vector
  signature_vector <- results_df$log2FoldChange
  names(signature_vector) <- toupper(results_df$gene_symbol)
  
  signature_vector <- signature_vector[!(duplicated(names(signature_vector)) | duplicated(names(signature_vector),fromLast=TRUE))]

  #We'll transform the vector into a simple -1/+1 vector
  signature_vector[which(signature_vector<0)] <- -1
  signature_vector[which(signature_vector>0)] <- 1
  
  #Now we create the signature, append it to the other signatures and then move on
  sigData <- na.omit(signature_vector)
  #We will now ask if the vector is longer than 20 - if not, there isnt much to build a signture from
  #We are now creating a signature, which is stored in the mySignatues bector
  if (length(signature_vector) > 20) {
    sig_tissue_aginggenes <- createGeneSignature(name = paste('Unique_', tissue,'_', 'UpDownReg_', 'AgingSignature', sep=''), sigData = sigData)
    mySignatures_CAS_AgingDEGs_perRegion <- c(mySignatures_CAS_AgingDEGs_perRegion, sig_tissue_aginggenes)
  }
  sig_tissue_aginggenes <- NULL
  sigData <- NULL
  
}

#We have now the signature library of common and region-specific DEGs assembled and will save this for calculatino of scores
save(mySignatures_CAS_AgingDEGs_perRegion, file='Signature_repository/AgingCohort_BulkDEG_signatures.bin')
