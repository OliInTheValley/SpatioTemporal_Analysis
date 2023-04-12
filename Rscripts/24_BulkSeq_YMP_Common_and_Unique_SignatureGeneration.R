
###Here, we want to inspect which DEGs are shared between tissues/regions and which are unique in response to a DR
#We load the DEseq2 object containing the bulk seq results for the aDR mice
load('R_objects/dds_BulkSeq_Rejuvenation_YMP.bin')


#First we assemble in which tissue which gene is passing the definition of an rejuvenation-related DEG (passing sinificance threshold)
#We set the significance threshold
padj_cutoff <- 0.05
#And the comparisons we want to use to define the DEGs
comparisons_to_analyze <-  c('YMP_vs_PBS')
#We set the number comparisons a gene needs to pass statistical testing - since we only have one comparison, we set this to 1
min_Difftimepoints <- 1
#We focus only on the results coming from analysing female YMP mice 
results_list_to_analyze <- results_list_CA2_YMP$Male  
#We will need to establish a list that we'll add the extracted information; here we'll add a table for each tissue where we summarise for each gene informatino gathered across the comparison YMP vs PBS
CA2_YMPDEGS_joined_list <- list()
#We set up a dummy dataframe onto which we will append results
comb_data_frame <- data.frame('DG'=1, row.names = 'SomeGene' )
#We exclude the hippocampus (posterior) from this analysis so that we won't be biased in favor or against hippocampus (as we have two tissue punches from that region)
#For this purpose we exclude it using the everything_but function
for (tissue in everything_but(names(results_list_to_analyze), c('hi2'))) {
  print(tissue)
  #In order for the later plotting to work, we will need tables with the exact same genes (i.e. all quanitfied ones) with some information which one was an rejuvenation-related DEG
  #We will iterate over the rejuvenation comparison (YMP vs PBS) and assemble the resultstables and later aggregate the information
  #We will need a dummy dataframe to append to. Therefore we simply take the first lane of the comparison YMP vs PBS and build basic dataframe and append
  #Only take genes with a padj beneath the cutoff to ensure no Deseq2 "NA" data is introduced 
  comparison_df <- head(subset(as.data.frame(results_list_to_analyze[[tissue]][[comparisons_to_analyze[1]]]$resall), padj <= 1), n=1)
  #We need a column in which we write what comparison this gene has been coming from. For this first line, we use a dummy value "init"
  comparison_df$tested_comparison <- 'init'
  for (comparison in comparisons_to_analyze) {
    #Well get the resultstable of that particular comparison and append
    comparison_df_temp <- as.data.frame(as.data.frame(results_list_to_analyze[[tissue]][[comparison]]$resall))
    comparison_df_temp$tested_comparison <- comparison
    comparison_df <- rbind(comparison_df, comparison_df_temp)
  }
  #First we remove the initial row
  comparison_df <- comparison_df %>% filter(tested_comparison != 'init')
  #Then we use dplyr to aggregate the information regarding in how many comparisons a gene has been crossing the significance threshould ('occured')
  #We'll also gather ifnormation about the  log2FC , as well as the baseMean. 
  temp_df <- comparison_df %>% 
    dplyr::group_by(gene_symbol) %>% 
    dplyr::summarise(occurance = sum(padj <= padj_cutoff), log2FoldChange=median(log2FoldChange), padj=min(padj) ,baseMean=median(baseMean))
  #Assign if a gene is an rejuvenation-related DEG if it passed the minimally requied number of comparisons the significance cutoff - in this case it's simply a 1
  temp_df <-  temp_df %>% mutate(RejuvenationRelated_DEG=ifelse(is.na(padj), 0, ifelse(occurance >=  min_Difftimepoints, ifelse(padj < padj_cutoff, 1, 0), 0) ))   
  #Likewise, we'll use the log2FC information to define if the gene is up/down-regualted or not affected (nm)
  temp_df$direction <- unlist(lapply(temp_df$log2FoldChange, function(x) ifelse(is.na(x), 'nm', ifelse(x > 0, 'up', 'down'))))
  temp_df <- as.data.frame(temp_df)
  #Well set the genes as row.names to fascilitate the later reformatting
  row.names(temp_df) <- temp_df$gene_symbol
  #We select and carry forward the columns that are critical to us
  temp_df <- temp_df[c('RejuvenationRelated_DEG', 'direction', 'occurance', 'gene_symbol', 'log2FoldChange', 'padj', 'baseMean')]
  #The result should be a dataframe that has always the same length (21,076 genes), ie all quantified genes. This is essential to be the case across all tissues
  CA2_YMPDEGS_joined_list[[tissue]] <- temp_df
  
}
#we use cbind.data.frame to collapse all dataframes from each tissue together
plot_df <- cbind.data.frame(CA2_YMPDEGS_joined_list)

#Now we'll extract all the columns that contain the ifnromation about wether the gene is an age-related DEG
#The reuslt will be a matrix that contains the tissues as columns and the rows corresponding to a gene, with values of 0/1 indicating if a gene is a DEG
plot_df_cond <- plot_df[,grepl('RejuvenationRelated_DEG',colnames(plot_df))]
#We'll have to reformat a bit
colnames(plot_df_cond) <- gsub('\\.RejuvenationRelated_DEG', '', colnames(plot_df_cond))
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
#First we'll subset this dataframe to genes that are at in at least 1 tissue an rejuvenation-related DEG
plot_df_cond_min <- plot_df_cond[which(rowSums(plot_df_cond[tissues_to_plot])>0),]
plot_df_cond_min <- as.data.frame(plot_df_cond_min)
#This is a very useful representation of all the DEG results so we'll label it properly and store it as an R object
CA2_YMPDEGS_genetable_min <- plot_df_cond_min

#The joind_list that contains for each tissue a single table that summarises the results from across all rejuvenation-realted comparisons will be useful in the future. We'll save it together  
#With the DEG_genetable_min

save(CA2_YMPDEGS_genetable_min, CA2_YMPDEGS_joined_list, file='R_objects/BulkSeq_RejuvenationYMPDEGS_genetable.bin')

#There should be only rows left where at least one tissue exhibits a 1 (ie rejuvenation-related DEG)
#We create a dataframe with boolean values
plot_df_cond_min_boolean <- CA2_YMPDEGS_genetable_min[tissues_to_plot] == 1

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
gene_count <- sort(rowSums(CA2_YMPDEGS_genetable_min), decreasing = T)
#We turn this into a dataframe and plot this as a bargarph
gene_count_df <- data.frame(count=gene_count, row.names = names(gene_count))
#Now we plot this and generate Figure 7G
ggplot(data=gene_count_df, aes(x=count)) + geom_bar() + scale_x_continuous(limits = c(0,15), expand = c(0,0), breaks = seq(1,15,1))

### Intersetingly, it seems that there's not much of a shared expression across regions. YMP appears to have more region-specific effects
#Something that's partiuclalry notable is strong differential expression found specifically in the SVZ. 
#When we check out some of the genes, it's notable that there are several neurodevelopmental genes induces, as well as genes related to neural stem cell differentiation
#This makes us wonder if there's some specific gene signature in the SVZ that we could map onto single cell data. 
#Let's create a couple of signatures of some of the most affected regions (SVZ, corpus callosum, caudate putamen) of both up/down-regualted genes, as well as signatures capturing only up or down-regulated genes


mySignatures_YMP_DEG <- c()

for (tissue in c('svz', 'cc', 'cp')) {
  #As always with VISION we need to supply which direcitno the gene is going (up/down) as demarcated by a +1/-1, respectively
  #We subset for the respective tissu
  results_df <- as.data.frame(CA2_YMPDEGS_joined_list[[tissue]])
  #Only consider genes that have been defined as 'RejuvenationRelated_DEG'
  results_df <- results_df %>% dplyr::filter(RejuvenationRelated_DEG == 1)
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
    sig_tissue_rejuvinationgenes <- createGeneSignature(name = paste(tissue,'_', 'UpDownReg_', 'RejuvenationSignature', sep=''), sigData = sigData)
    mySignatures_YMP_DEG <- c(mySignatures_YMP_DEG, sig_tissue_rejuvinationgenes)
  }
  sig_tissue_rejuvinationgenes <- NULL
  sigData <- NULL
  
  #Now we'll only build a signature out of the up-regualted genes
  signature_vector_upregulated <- subset(signature_vector, signature_vector >0)
  signature_vector_downregulated <- subset(signature_vector, signature_vector <0)
  
  #Now we create the up-regulated signature, append it to the other signatures and then move on
  sigData <- na.omit(signature_vector_upregulated)
  #We will now ask if the vector is longer than 20 - if not, there isnt much to build a signture from
  #We are now creating a signature, which is stored in the mySignatues bector
  if (length(signature_vector_upregulated) > 20) {
    sig_tissue_rejuvinationgenes <- createGeneSignature(name = paste(tissue,'_', 'UpReg_', 'RejuvenationSignature', sep=''), sigData = sigData)
    mySignatures_YMP_DEG <- c(mySignatures_YMP_DEG, sig_tissue_rejuvinationgenes)
  }
  sig_tissue_rejuvinationgenes <- NULL
  sigData <- NULL
  
  
  #Now we create the down-regulated signature, append it to the other signatures and then move on
  #VISION, however, requires that there can only be a mixed-sign or all-positive-sign vector. So we have to turn the signs on the down-regulated vector
  signature_vector_downregulated <- signature_vector_downregulated*-1
  sigData <- na.omit(signature_vector_downregulated)
  #We will now ask if the vector is longer than 20 - if not, there isnt much to build a signture from
  #We are now creating a signature, which is stored in the mySignatues bector
  if (length(signature_vector_downregulated) > 20) {
    sig_tissue_rejuvinationgenes <- createGeneSignature(name = paste(tissue,'_', 'DownReg_', 'RejuvenationSignature', sep=''), sigData = sigData)
    mySignatures_YMP_DEG <- c(mySignatures_YMP_DEG, sig_tissue_rejuvinationgenes)
  }
  sig_tissue_rejuvinationgenes <- NULL
  sigData <- NULL
  
}

#We will save this for calculatino of scores
save(mySignatures_YMP_DEG, file='Signature_repository/RejuvenationCohort_YMP_signatures.bin')
