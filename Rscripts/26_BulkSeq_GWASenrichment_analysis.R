##In this final analysis, we will perform a type of functional enrichment analysis fpr DEGs found in each region
#Instead of GO-terms, we will use gene lists of human GWAS hits for three neurodegenerative disorders: Alzheimer's disease (AD), Prakinson's disease (PD) and Multiple Sclerosis.
#AD and PD are well-known age-related diseases of the brain and MS is a well-accepted white matter disease, which we consider interesting in the context of our study due to the accelrated aging trajctories in white matter-rich areas

#We will perform the enrichment analysis for all three diseases and all regions except of the posterior hippocampus (we will use the anterior hippocampus as equivaelent), the olfacotry bulb and the entorhinal cortex - the latter exhibited too few age-related DEGs for a reasonable anlaysis to work

#We load our Deseq2 results
load('R_objects/dds_BulkSeq_Aging.bin')

#We will use the list of GWAS hits that was utilized in a previously publishe study: Yang & Kern et al. 2021, doi: 10.1038/s41586-021-03710-0
GWAS_trait_human_disease_genes_curated_mouse_mapped_long_format <- read_delim("input_data/GWAS_enrichment/GWAS_trait_human_disease_genes_curated_mouse_mapped_long_format.txt", 
                                                                              delim = "\t", escape_double = FALSE, 
                                                                              trim_ws = TRUE)
#We filter out the three regions mentioned above
tissues_to_analyze <- everything_but(names(results_list_CA1$Both), c( 'hi2', 'olf', 'ent'))

#We will need qutie some information to perform the enrichment analysis
#We thus prepare respective dataframes to hold that
data_frame_total <- data.frame(init=rep(1, length(tissues_to_analyze)), row.names = tissues_to_analyze)
data_frame_fisher_pval <- data.frame(init=rep(1, length(tissues_to_analyze)), row.names = tissues_to_analyze)
data_frame_fisher_oddsRatio <- data.frame(init=rep(1, length(tissues_to_analyze)), row.names = tissues_to_analyze)
data_frame_fisher_stars <- data.frame(init=rep('*', length(tissues_to_analyze)), row.names = tissues_to_analyze)

#this is a lis where we'll store which GWAS hit of which disease is a DEG in a given tissue
common_geneList <- list()


data_frame_collection <- data.frame(category='INIT', tissue='INIT', DEGs=c('Gene1, Gene2'), stringsAsFactors = F)
i <- 1
#We will rotate over each disease indepently and perform the same type of analysis steps
for (category in unique(GWAS_trait_human_disease_genes_curated_mouse_mapped_long_format$Disease)) {
  print(category)
  common_geneList[[category]] <- list()
  #First, we extract the GWAS hits associated with the current disease
  gene_list_to_check <- GWAS_trait_human_disease_genes_curated_mouse_mapped_long_format %>% filter(Disease == category) %>% pull(Gene)
  #Now we'll go over each tissue and see a) which of these GWAS hits are expressed - as defined by Deseq2's independent filtering step - and b) which GWAS hits are differentially expressed
  for (tissue in tissues_to_analyze) {
    #As always, we'll get our list of age-related DEGs, which are classified as crossing the padj threshold of 0.05 in at least two of the comparisons between 3 months and any following age group
    DEGs_of_tissue <- (unlist((lapply(results_list_CA1$Both[[tissue]][c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3', '28_vs_3' )], function (x) (row.names(subset(x[["resall"]], padj < 0.05)))))))
    #If a gene occurs in more than two comparisons, we define it as DEG
    DEGs_of_tissue <- names(which(table(DEGs_of_tissue) >= 2))
    #As for any functional enrichment analysis, we'll require a list of 'expressed' genes as background, which we define according to Deseq2's definiton (a gene is expressed when it has a non-NA value in the padj column)
    #We'll use the 3 vs 21 months comparison for that to span a range of ages in our dataset
    Expressed_in_tissue <- results_list_CA1$Both[[tissue]]$`21_vs_3`$resOrdered %>% filter(!is.na(padj)) %>% pull(gene_symbol)

    #We'll subset our GWAS hits for those that are expressed in the region/tissue
    gene_list_to_check_inTissue <- intersect(gene_list_to_check, Expressed_in_tissue)
    
    #Now we count the expressed GWAS hits that are age-related DEGS
    common_genes <- length(intersect(gene_list_to_check_inTissue, DEGs_of_tissue))
    #We'll save that information in the common_geneList
    common_geneList[[category]][[tissue]] <- intersect(gene_list_to_check_inTissue, DEGs_of_tissue)
    #We will need to know which GWAS hits are NOT DEGs
    uncommon_diseasegenes <- length(setdiff(gene_list_to_check_inTissue, DEGs_of_tissue))
    #We will need to know which DEGs hits are NOT GWAS hits
    uncommon_DEGs <- length(setdiff(DEGs_of_tissue, gene_list_to_check_inTissue))
    #From these four 'intredients' (GWAS/DEG overlap, DEG-only, GWAS-only, background expressed genes) we'll build the hypergeometric test for the enrichement test
    test_matrix <- matrix(c(common_genes, uncommon_diseasegenes, uncommon_DEGs, length(Expressed_in_tissue)-uncommon_diseasegenes - uncommon_DEGs - common_genes), ncol = 2)
    #Now we'll run Fisher's exact test and extract the p value.  As we're only interested in enrichemnts, we perform this as a one-sided test
    fisher_testresult <- fisher.test(test_matrix, alternative = 'greater')
    p_hypergeometric <- fisher_testresult$p.value
    #We'll store the pvalue and odds ratio estimate
    data_frame_fisher_pval[tissue,category] <- p_hypergeometric
    data_frame_fisher_oddsRatio[tissue,category] <- fisher_testresult$estimate
    #We'll also require information about the number of GWAS-DEGs as well as the 
    data_frame_total[tissue,category] <- length(intersect(gene_list_to_check_inTissue, DEGs_of_tissue))
    i <- i + 1
    data_frame_collection[i, 'category'] <- category
    data_frame_collection[i, 'tissue'] <- tissue
    data_frame_collection[i, 'DEGs'] <- paste(intersect(gene_list_to_check_inTissue, DEGs_of_tissue), collapse = ', ')
    
  }
  
}

#
data_frame_total$init <- NULL
data_frame_fisher_pval$init <- NULL
data_frame_fisher_oddsRatio$init <- NULL
data_frame_fisher_stars$init <- NULL

#We perform multiple testing correction and convert the pvalues into astericks to inspect in which region we found a significant enrichment
data_frame_fisher_pval_t <- as.data.frame(t(data_frame_fisher_pval))
data_frame_fisher_stars <- as.data.frame(data_frame_fisher_pval_t)
for (tissue in colnames(data_frame_fisher_pval_t)) {
  data_frame_fisher_pval_t[,tissue] <- p.adjust(data_frame_fisher_pval_t[,tissue], method = 'BH')
  data_frame_fisher_stars[,tissue] <- apply(data_frame_fisher_pval_t[tissue], 1, stars.return)
}


#Having idenitfied the regions with significant enrichment, we'll create a so-called arc plot visually reprents the region's individual DEG enrichment as well as their potential overlap via arcs
###We will exhibit this for MS as an example, and this can be repeated for AD and PD, respectively
#First, we'll have to reshpae the stats dataframe and the GWAS-DEG count table
plot_df_stats <- melt(as.matrix(data_frame_fisher_pval_t))
colnames(plot_df_stats) <- c('Disease','Region', 'pval')
plot_df_stats$region_disease <- paste(plot_df_stats$Region, plot_df_stats$Disease, sep = '_')

plot_df_count <- melt(t(data_frame_total))
colnames(plot_df_count) <- c('Disease','Region', 'count')
plot_df_count$region_disease <- paste(plot_df_count$Region, plot_df_count$Disease, sep = '_')
#Same thing for the odds ratio estimates. This will form the plot backbone, and we'll attach ohter information from the stats and count data table over
plot_df <- melt(as.matrix(t(data_frame_fisher_oddsRatio)))
colnames(plot_df) <- c('Disease', 'Region','enrichment')
plot_df$region_disease <- paste(plot_df$Region, plot_df$Disease, sep = '_')
plot_df$pval <- as.numeric(plyr::mapvalues(plot_df$region_disease, from = plot_df_stats$region_disease, to = plot_df_stats$pval))
plot_df$count <- as.numeric(plyr::mapvalues(plot_df$region_disease, from = plot_df_count$region_disease, to = plot_df_count$count))

#Now we'll start the plotting for MS 
#We will need to have the CellPlot package installed and loaded for that
category <- 'Multiple sclerosis'
#We subset for the disease we want to plot - in this case MS. We also will plot the Arc plot for all regions even if they do not exhibit a significant enrichment
#This way we'll have plots of the same shape and form across diseases. We'll highlight the ones with significant enrichment in the resulting plot later 
plot_df_Disease <- plot_df %>% filter(Disease == category) 
#We'll need the respective GWAS genes 
gene_list_to_check <- GWAS_trait_human_disease_genes_curated_mouse_mapped_long_format %>% filter(Disease == category) %>% pull(Gene)
#Arc plot is designed to work with GO terms, so 
colnames(plot_df_Disease)[2] <- 'Region'
colnames(plot_df_Disease)[3] <- 'OddsRatio'

#we will need to know which GWAS/DEGs are up/down regulated to build the resulting plot
genes_up = list()
for (current_region in plot_df_Disease$Region){
  #We have to - once again - identify the GWAS/DEGs. So we'll get the DEGs, the expressed genes, and GWAS hits
  DEGs_of_tissue <- (unlist((lapply(results_list_CA1$Both[[current_region]][c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3', '28_vs_3' )], function (x) (row.names(subset(x[["resall"]], padj < 0.05)))))))
  DEGs_of_tissue <- names(which(table(DEGs_of_tissue) >= 2))
  Expressed_in_tissue <- results_list_CA1$Both[[current_region]]$`21_vs_3`$resOrdered %>% filter(!is.na(padj)) %>% pull(gene_symbol)

  gene_list_to_check_inTissue <- intersect(gene_list_to_check, Expressed_in_tissue)
  #This gives us the expressed GWAS hits that are DEGs in the current tissue
  common_genes <- intersect(gene_list_to_check_inTissue, DEGs_of_tissue)
  #We now need the direction of regulation which we get from the 21 vs 3 months comparison - we only take the genes that are up-regulated
  expression_temp <- results_list_CA1$Both[[current_region]]$`21_vs_3`$resOrdered %>% filter(gene_symbol %in% common_genes) %>% filter(log2FoldChange > 0) %>% pull(gene_symbol)
  genes_up = c(genes_up, list(expression_temp))
}
names(genes_up) <- plot_df_Disease$Region

#We repeat all this for the down-regulated genes
genes_down = list()
for (current_region in plot_df_Disease$Region){
  
  DEGs_of_tissue <- (unlist((lapply(results_list_CA1$Both[[current_region]][c('12_vs_3', '15_vs_3', '18_vs_3', '21_vs_3', '26_vs_3', '28_vs_3' )], function (x) (row.names(subset(x[["resall"]], padj < 0.05)))))))
  DEGs_of_tissue <- names(which(table(DEGs_of_tissue) >= 2))
  Expressed_in_tissue <- results_list_CA1$Both[[current_region]]$`21_vs_3`$resOrdered %>% filter(!is.na(padj)) %>% pull(gene_symbol)

  gene_list_to_check_inTissue <- intersect(gene_list_to_check, Expressed_in_tissue)
  
  common_genes <- intersect(gene_list_to_check_inTissue, DEGs_of_tissue)
  expression_temp<- results_list_CA1$Both[[current_region]]$`21_vs_3`$resOrdered %>% filter(gene_symbol %in% common_genes) %>% filter(log2FoldChange < 0) %>% pull(gene_symbol)
  genes_down = c(genes_down, list(expression_temp))
}
names(genes_down) <- plot_df_Disease$Region

#Now we have all the ingredients for building the arcplot 
#This creates the basis for Figure panels 8E-G
arc.plot(x = setNames(plot_df_Disease$OddsRatio, plot_df_Disease$Region), 
         up.list = genes_up, 
         down.list = genes_down, 
         x.mar = c(0.8, 1), t=0.25, tc = 'black', fixed.scale = 1, x.bound = 2)

#The arcplot produced from this still requires some manual filtering, espcially removing of arcs that connects significantly-enriched from non-significantly enriched regions
#A simplified version of this plot can be created if filtering the plot_df for significantly enriched regions only


