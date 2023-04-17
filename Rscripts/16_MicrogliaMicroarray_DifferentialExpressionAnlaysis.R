#Here, we will perfom a meta analysis of microarray data of bulk sorted microglia from multiple regions
#The data was published by Garbert et al., 2016 DOI:10.1038/nn.4222 
################################################################
#   Here we peform the differential expression analysis with limma
#   load series and platform data from GEO

gset <- getGEO("GSE62420", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL11180", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# then we group membership for all samples
gsms <- "AAAABBBBCCCCDDDDXXXXXXXXEEEEFFFFGGGGHHHHIIIIJJJJKKKKLLLL"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Ceb_04M","Cor_04M","Hip_04M","Str_04M",
                       "Ceb_12M","Cor_12M","Hip_12M","Str_12M",
                       "Ceb_22M","Cor_22M","Hip_22M","Str_22M"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=T)
# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
# we are mainly interested in age-related expression effects within each region - that is, we won't explore inter-regino differences
cts <- c("Ceb_12M-Ceb_04M",  "Str_12M-Str_04M", "Hip_12M-Hip_04M", "Cor_12M-Cor_04M",  
         "Ceb_22M-Ceb_04M",  "Str_22M-Str_04M", "Hip_22M-Hip_04M", "Cor_22M-Cor_04M",
         "Ceb_22M-Ceb_12M",  "Str_22M-Str_12M", "Hip_22M-Hip_12M", "Cor_22M-Cor_12M",
         "Ceb_22M-Cor_22M",  "Str_22M-Cor_22M", "Ceb_22M-Hip_22M", "Str_22M-Hip_22M")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results for age comparison 12 vs 3m
vennDiagram(dT[,1:4], circle.col=palette())
# Venn diagram of results for age comparison 22 vs 3m
vennDiagram(dT[,5:8], circle.col=palette())
# we can apprecaite that there are more differentially expressed genes in the striatum and cerebellum in particular at mid-age

# Bargraph of DEGs by up/downregulation

plot_df <- data.frame(n_genes= c(table(dT[,1])[1], table(dT[,2])[1], table(dT[,3])[1], table(dT[,4])[1], 
                                 table(dT[,5])[1], table(dT[,6])[1], table(dT[,7])[1], table(dT[,8])[1],
                                 table(dT[,1])[3], table(dT[,2])[3], table(dT[,3])[3], table(dT[,4])[3], 
                                 table(dT[,5])[3], table(dT[,6])[3], table(dT[,7])[3], table(dT[,8])[3]),
                      region= c('ceb', 'str', 'hip', 'cor',
                                'ceb', 'str', 'hip', 'cor',
                                'ceb', 'str', 'hip', 'cor',
                                'ceb', 'str', 'hip', 'cor'),
                      direction= c('up', 'up', 'up', 'up',
                                   'up', 'up', 'up', 'up',
                                   'down', 'down', 'down', 'down',
                                   'down', 'down', 'down', 'down'),
                      comparison= c('M12', 'M12', 'M12', 'M12',
                                    'M22', 'M22', 'M22', 'M22',
                                    'M12', 'M12', 'M12', 'M12',
                                    'M22', 'M22', 'M22', 'M22')
)

# order regions vector
plot_df$region <- factor(plot_df$region, levels = c('ceb', 'str', 'hip', 'cor'), ordered = T)
plot_df$direction <- factor(plot_df$direction,levels = c('up', 'down'), ordered = T)

#This creates Figure S12D
myplot <- ggplot(plot_df , aes(x=region, y=n_genes, fill=direction), color=NA) + geom_bar(stat = "identity") + 
   xlab("") +   
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1)) + facet_grid(.~comparison, scales = 'free_y') 
myplot


# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# Now we'll plot the log2FC distribution (old vs young) for each region for the CAS genes
# We have to focus here on the ones that become up-regulated with age, as we do not have enough CAS genes that become down-regulated to perform a statistically solid analysis for them
# So we load the CAS genes and the results tables from the corpus callosum to extract the direction of thre regulation

load('R_objects/CA1_CASgenes.bin')
load('R_objects/dds_BulkSeq_Aging.bin')

genes_to_analyze <- results_list_CA1$Both$cc$`26_vs_3`$resOrdered[CA1_CASgenes,] %>% filter(log2FoldChange >0) %>% pull(gene_symbol)

#We also use the infromation in fit2 to extract all probe IDs that match to the CAS genes - we do acknowledge that individual CAS genes will have no probe
annotatino_df <- fit2$genes
#Now we extract the respective probes
probes_to_analyze <- annotatino_df %>% filter(Gene.symbol %in% genes_to_analyze)  %>%  pull(ID)

#We can further reduce the probe number by filerting for probes that exhibit a significant change in at least one of the 8 age-related comparisons (i.e. 12 vs 4 months or 22 vs 4 months)
#This is optional, as the results stay the same. In the manuscript, we have have performed this step
probes_to_analyze <- intersect(probes_to_analyze, row.names(as.data.frame(dT)[rowSums(abs(as.data.frame(dT)[,c(1:8)])) >= 1,]))

log2table <- as.data.frame(fit2$coefficients)

log2table$probeName <- row.names(log2table)
#Now we can plot the distribution of log2FCs with age in each region
#This creates the Figure S12 H
(myplot <- ggplot(data=log2table[probes_to_analyze,], aes(x='None',  y=(log2table[probes_to_analyze,comparison])) ) + 
    geom_violin(aes(x='Ceb', y=log2table[probes_to_analyze,'Ceb_12M-Ceb_04M']), color='red',draw_quantiles = c(0.5)) +
    geom_violin(aes(x='Str', y=log2table[probes_to_analyze,'Str_12M-Str_04M']), color='purple',draw_quantiles = c(0.5)) +
    geom_violin(aes(x='Hip', y=log2table[probes_to_analyze,'Hip_12M-Hip_04M']), color='black',draw_quantiles = c(0.5)) +
    geom_violin(aes(x='Cor', y=log2table[probes_to_analyze,'Cor_12M-Cor_04M']), color='green',draw_quantiles = c(0.5)) +
    geom_jitter(aes(x='Ceb', y=log2table[probes_to_analyze,'Ceb_12M-Ceb_04M']), color='black') +
    geom_jitter(aes(x='Str', y=log2table[probes_to_analyze,'Str_12M-Str_04M']), color='black') +
    geom_jitter(aes(x='Hip', y=log2table[probes_to_analyze,'Hip_12M-Hip_04M']), color='black') +
    geom_jitter(aes(x='Cor', y=log2table[probes_to_analyze,'Cor_12M-Cor_04M']), color='black') +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-1,1), expand = c(0,0)) 
)

#And we can test for significant differences between regions by setting up a dataframe where we can perform paired, pairwise wilcoxon tests
table_anova_input <- log2table[probes_to_analyze,c('Ceb_12M-Ceb_04M', 'Str_12M-Str_04M', 'Hip_12M-Hip_04M', 'Cor_12M-Cor_04M')]
colnames(table_anova_input) <- c('ceb', 'str', 'hip', 'cor')
table_anova_input <- melt(as.matrix(table_anova_input))

colnames(table_anova_input) <- c('probe', 'region', 'log2FC')
pwc <- table_anova_input %>%
  pairwise_wilcox_test(
    log2FC ~ region, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc


#We can repeat this analysis for the comparison 22 vs 4 months, which yields the other half of the Figure S12H

(myplot <- ggplot(data=log2table[probes_to_analyze,], aes(x='test', y=(log2table[probes_to_analyze,comparison])) ) + 
    geom_violin(aes(x='Ceb', y=log2table[probes_to_analyze,'Ceb_22M-Ceb_04M']), color='red',draw_quantiles = c(0.5)) +
    geom_violin(aes(x='Str', y=log2table[probes_to_analyze,'Str_22M-Str_04M']), color='purple',draw_quantiles = c(0.5)) +
    geom_violin(aes(x='Hip', y=log2table[probes_to_analyze,'Hip_22M-Hip_04M']), color='black',draw_quantiles = c(0.5)) +
    geom_violin(aes(x='Cor', y=log2table[probes_to_analyze,'Cor_22M-Cor_04M']), color='green',draw_quantiles = c(0.5)) +
    geom_jitter(aes(x='Ceb', y=log2table[probes_to_analyze,'Ceb_22M-Ceb_04M']), color='black') +
    geom_jitter(aes(x='Str', y=log2table[probes_to_analyze,'Str_22M-Str_04M']), color='black') +
    geom_jitter(aes(x='Hip', y=log2table[probes_to_analyze,'Hip_22M-Hip_04M']), color='black') +
    geom_jitter(aes(x='Cor', y=log2table[probes_to_analyze,'Cor_22M-Cor_04M']), color='black') +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-1,2.5), expand = c(0,0), oob=squish) 
)

table_anova_input <- log2table[probes_to_analyze,c('Ceb_22M-Ceb_04M', 'Str_22M-Str_04M', 'Hip_22M-Hip_04M', 'Cor_22M-Cor_04M')]
colnames(table_anova_input) <- c('ceb', 'str', 'hip', 'cor')
table_anova_input <- melt(as.matrix(table_anova_input))

colnames(table_anova_input) <- c('probe', 'region', 'log2FC')
pwc <- table_anova_input %>%
  pairwise_wilcox_test(
    log2FC ~ region, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc


