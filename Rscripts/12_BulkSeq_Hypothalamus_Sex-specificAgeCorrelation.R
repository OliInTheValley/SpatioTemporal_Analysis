
rm(list=ls())

##Here, we perform will correlate expression for each gene in the hypothalamus with age for males and females separately and assess to what degree the expression data overlaps between the two sexes
#Importantly, we'll have to do this for the time window of 3 to 21 months, as we only have data from both sexes for that time frame
#First, we load our Deseq2 object as it contains all the counts and metadata
load('R_objects/dds_BulkSeq_Aging.bin')
load('R_objects/CA1_CASgenes.bin')

#We'll setup a results list in which we'll store the output of the following analysis
BulkSeq_Agecorrelation_results_list <- list()
#We'll set the tissue to hypothalamus
tissue <- 'hy'
#We extract the respective Deseq2 object
dds_hy <- dds_CA1_list[[tissue]]

#First we subset the deseq2 object to ages 3 to 21 months
dds_hy <- dds_hy[, dds_hy$age %in% c(3,12,15,18,21) ]
dds_hy$age <- droplevels(dds_hy$age)
#Now we create a subset of only females and males
dds_hy_male  <- dds_hy[, dds_hy$sex == 'Male']
dds_hy_female  <- dds_hy[, dds_hy$sex == 'Female']

#We'll drop unused factor levels from the sex column 
dds_hy_male$sex <- droplevels(dds_hy_male$sex)
dds_hy_female$sex <- droplevels(dds_hy_female$sex)

#We'll alter the design matrix and then run Deseq2
design(dds_hy_male) <- ~age
dds_hy_male <- DESeq(dds_hy_male, fitType = 'local')
design(dds_hy_female) <- ~age
dds_hy_female <- DESeq(dds_hy_female, fitType = 'local')


#To make sure that we are not correlating flat noise with flat noise, we will focus on genes that
#we consider as 'expressed', which means we'll discard all genes that Deseq2 has flagged in their independent filtering by assigning a na in the padj column
#We'll use the results list from 3 to 21 months to make sure we capture genes only expressed at old age
results_male <- results(dds_hy_male, contrast = c('age', '21', '3'), cooksCutoff = T) #Set contrast and retrieve results
results_male$gene_symbol <- row.names(results_male) #get the rownames as a separate column in the results report
resOrdered_male <- as.data.frame(results_male[order(results_male$pvalue),]) #create a simple dataframe ordered by the padj
#And we extract the set of expressed genes in males
ExpresedGenes_male <- resOrdered_male %>% filter(!is.na(padj)) %>% pull(gene_symbol)

#We'll repeat that for females
results_female <- results(dds_hy_female, contrast = c('age', '21', '3'), cooksCutoff = T) #Set contrast and retrieve results
results_female$gene_symbol <- row.names(results_female) #get the rownames as a separate column in the results report
resOrdered_female <- as.data.frame(results_female[order(results_female$pvalue),]) #create a simple dataframe ordered by the padj
#And we extract the set of expressed genes in females
ExpresedGenes_female <- resOrdered_female %>% filter(!is.na(padj)) %>% pull(gene_symbol)

#Interestingly, there are a few more genes detected as 'expressed' in females than in males
length(ExpresedGenes_male)
length(ExpresedGenes_female)
#It is possible that there are few genes that become expressed with age in either sex that would never appear in the other sex. So to not bias 
#the analysis against one or the other sex, we'll merge these two sets of expressed genes and use them as input for the correlation analyses
ExpressedGenes_inEitherSex <- unique(c(ExpresedGenes_male, ExpresedGenes_female))

##Now we'll run the age-correlation analysis for each gene for the male samples
#We use the DESE2-normalized data, which we extract by getting the counts from the Deseq2 object
count_table <- counts(dds_hy_male, normalized=T)
count_table <- as.data.frame(count_table)
#We definie the numberic values we want to correlate against - in this case, the age of the mice
Score_to_compare <- as.numeric(as.character(as.data.frame(colData(dds_hy_male))[,'age']))

#Perform gene-vs-age score correlataion
#To do that we set up a dataframe that can holde the results
#We want to record both the spearman and pearson correlation corefficient, as well as the results from cor.test to see that we have statistical
#confidence in our correlation coefficents
correlation_frame <- data_frame(gene='A', cor_spear=1, p_spear=1, cor_pears=1, p_pears=1)

#We set up a counter so that we can have a progressbar for our analysis
i <- 1
pb = txtProgressBar(min = 1, max = length(ExpressedGenes_inEitherSex), initial = 1, style = 3)
#We will now iterate over very expresed gene and perform the following steps
for (gene in ExpressedGenes_inEitherSex) {
  setTxtProgressBar(pb,i)
  #We set up a temporary dataframe that holds age and normalized counts 
  cor_frame <- data.frame(expr=t(count_table[gene,]), concentration= Score_to_compare)
  colnames(cor_frame) <- c('expr', 'concentration')
  #We wil only retain samples where the count is not 0
  cor_frame <- cor_frame %>% filter(expr > 0)
  #After filtering for non-expressing samples, we will only continue the analysis if we have at least 15 samples left - otherwise we end up
  #with very spurious correlations. This is a bit of a lower threshold than what we applied for the age-correlation analysis using both sexes,
  #though we have to consider that there are fewer samples to begin with as we only focus on one sex and ages 3 to 21 months
  if (nrow(cor_frame) >= 15) {
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
    
    #if the gene had fewer than 15 samples expressing it, we'll add a dummy row
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
#correlation_frame <- na.omit(correlation_frame)
correlation_frame$padj_spear <- p.adjust(correlation_frame$p_spear, method = 'BH')
correlation_frame$padj_pear <- p.adjust(correlation_frame$p_pears, method = 'BH')
#Relabel
hy_correlation_fullTable_male <- correlation_frame
rm(correlation_frame)
#Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibhitng an absolute correlation value of 0.5)
print(paste('TopCorGenes: ', nrow(hy_correlation_fullTable_male %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
nrow(hy_correlation_fullTable_male %>% filter(padj_spear < 0.05))
#Filter and store data for high correlates 
hy_correlation_fullTable_topCor_male  <- (hy_correlation_fullTable_male %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))



##Now we'll run the age-correlation analysis for each gene for the female samples
#We use the DESE2-normalized data, which we extract by getting the counts from the Deseq2 object
count_table <- counts(dds_hy_female, normalized=T)
count_table <- as.data.frame(count_table)

#We definie the numberic values we want to correlate against - in this case, the age of the mice
Score_to_compare <- as.numeric(as.character(as.data.frame(colData(dds_hy_female))[,'age']))

#Perform gene-vs-age score correlataion
#To do that we set up a dataframe that can holde the results
#We want to record both the spearman and pearson correlation corefficient, as well as the results from cor.test to see that we have statistical
#confidence in our correlation coefficents
correlation_frame <- data_frame(gene='A', cor_spear=1, p_spear=1, cor_pears=1, p_pears=1)

#We set up a counter so that we can have a progressbar for our analysis
i <- 1
pb = txtProgressBar(min = 1, max = length(ExpressedGenes_inEitherSex), initial = 1, style = 3)
#We will now iterate over very expresed gene and perform the following steps
for (gene in ExpressedGenes_inEitherSex) {
  setTxtProgressBar(pb,i)
  #We set up a temporary dataframe that holds age and normalized counts 
  cor_frame <- data.frame(expr=t(count_table[gene,]), concentration= Score_to_compare)
  colnames(cor_frame) <- c('expr', 'concentration')
  #We wil only retain samples where the count is not 0
  cor_frame <- cor_frame %>% filter(expr > 0)
  #After filtering for non-expressing samples, we will only continue the analysis if we have at least 15 samples left - otherwise we end up
  #with very spurious correlations. This is a bit of a lower threshold than what we applied for the age-correlation analysis using both sexes,
  #though we have to consider that there are fewer samples to begin with as we only focus on one sex and ages 3 to 21 months
  if (nrow(cor_frame) >= 15) {
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
    
    #if the gene had fewer than 15 samples expressing it, we'll add a dummy row
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
hy_correlation_fullTable_female <- correlation_frame
rm(correlation_frame)
#Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibhitng an absolute correlation value of 0.5)
print(paste('TopCorGenes: ', nrow(hy_correlation_fullTable_female %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
nrow(hy_correlation_fullTable_female %>% filter(padj_spear < 0.05))
#Filter and store data for high correlates 
hy_correlation_fullTable_topCor_female  <- (hy_correlation_fullTable_female %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))


#We have correlation results for each of the sexes. We'll know combine these and inspect how these look like in comparison
comparison_list <- list(male=subset(results_male, padj < 0.1)$gene_symbol,
                        female=subset(results_female, padj < 0.1)$gene_symbol)

#Let's look at the direct overlap using vennerable's venn function
comparison_list <- list(male=hy_correlation_fullTable_topCor_male$gene,
                        female=hy_correlation_fullTable_topCor_female$gene)
plot(Venn(comparison_list))
#The direct intersect is not particularly strong but not necessarily unexpected in this tissue, given the significant difference in CAS scores, so we expexct the tissue
#to have differeing biological age in males/females
#Let's see what's in the intersectx
(commom_genes <- Reduce(intersect, comparison_list))
#This looks good. We see genes like C4b, Ctss, Gpr17, C1qb etc. so classcial CAS genes
#We will extract the genes that we found to correlate primarily in females or males, respectively
enriched_genes_female <- setdiff(hy_correlation_fullTable_topCor_female$gene, 
                               hy_correlation_fullTable_topCor_male$gene)

enriched_genes_male <- setdiff(hy_correlation_fullTable_topCor_male$gene,
                             hy_correlation_fullTable_topCor_female$gene)

#Now we merge the correlation data from each sex by the gene and mark which genes are commonly de-regulated, only detected in females or males.
gene_correlate_comparison_df <- merge.data.frame(hy_correlation_fullTable_female,
                                                 hy_correlation_fullTable_male, by = 'gene')

gene_correlate_comparison_df$topCorDEGrelate <- 'ns'
gene_correlate_comparison_df[gene_correlate_comparison_df$gene %in% enriched_genes_female, 'topCorDEGrelate'] <- 'female'
gene_correlate_comparison_df[gene_correlate_comparison_df$gene %in% enriched_genes_male, 'topCorDEGrelate'] <- 'male'
gene_correlate_comparison_df[gene_correlate_comparison_df$gene %in% commom_genes, 'topCorDEGrelate'] <- 'common'

#We'll also annotate any gene in the list that is part of the CAS set 
gene_correlate_comparison_df[gene_correlate_comparison_df$gene %in% CA1_CASgenes, 'CAS_gene'] <- 'no'
gene_correlate_comparison_df[gene_correlate_comparison_df$gene %in% CA1_CASgenes, 'CAS_gene'] <- 'yes'

#Now we'll plot a scatterplot of genes that show a significant (padj < 0.05 and absolute Spearman's coefficient ≥ absolute(0.5)) in either or both sex
#We'll drop all the genes flagged as non-sigfificant as this ends up associating noise with noise. 
#This creates panel S7C
(myplot <- ggplot(gene_correlate_comparison_df %>% filter(topCorDEGrelate != 'ns') , aes(x=cor_spear.x, y=cor_spear.y, fill=topCorDEGrelate, color=CAS_gene)) + 
    geom_point(shape=21)  
)

