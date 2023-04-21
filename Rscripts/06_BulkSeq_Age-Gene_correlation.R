#In this script, we will perform two different analyses to verify the results from the DEG analysis performed with DESeq2
#The first is correlating expression directly with age the other is running WGCNA and associating the resulting modules with changes with age
#Ideally, both analyses should reveal that certain regions (like the corpus callosum) exhibit many, accentuated shifts with age, while other regions, like the enthorinal cortex, exhibit very few effects
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





###Now we'll shift gears to the WGCNA approach
##########################################################################################################
###
###               WGCNA Co-expression Analysis  
###                   
###                Brain Aging Atlas
###
###            Author: Patricia Moran Losada. pmlosada@stanford.edu
###
##########################################################################################################
rm(list=ls())
dir.create('/WGCNA')
setwd("./WGCNA/")

options(stringsAsFactors = F)


##########################################################################################################
##              Input data
##########################################################################################################

#Load dataset from raw counttable

counts <- read_delim("input_data/BulkSeq_Aging_Counttable.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
counts <- as.data.frame(counts)
row.names(counts) <- counts[,1] #turn gene names into rownames of the table
counts[1] <- NULL

metadata <- read_delim("input_data/BulkSeq_Aging_metadata.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
metadata <- as.data.frame(BulkSeq_Aging_metadata)

##########################################################################################################
##             WGCNA per tissue
##########################################################################################################

acronym_labels <- c("cc", "cer","cor", "cp", "ent", "hi2","hi", "hy", "med", "olf", "plx",
                    "pon", "svz","th", "vis")

for(i in 1:length(acronym_labels)){
  metadata_tmp <- metadata[which(metadata$tissue %in% acronym_labels[i]),]
  counts_tmp <- counts[,which(colnames(counts) %in% metadata_tmp)] 
  
  table(rownames(metadata_tmp)==colnames(counts_tmp))
  counts_tmp <- counts_tmp[,rownames(metadata_tmp)]
  table(rownames(metadata_tmp)==colnames(counts_tmp))
  
  design=model.matrix(~age+sex,data = metadata_tmp)
  colnames(design)[1]="Intercept"
  dge.voom = voom(calcNormFactors(DGEList(counts = counts_tmp),method = 'TMM'), design, plot = T)
  datExpr = as.matrix(dge.voom$E) #normalized expression values on the log2 scale
  
  #check for excesive missing values
  gsg = goodSamplesGenes(datExpr, verbose = 3);
  gsg$allOK
  
  ### cluster
  sampleTree = hclust(dist(t(datExpr)), method = "average");
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 3,cex=2)
  
  
  ##network construction and module detection. Choose a set of soft-thresholding powers
  powers = c(seq(1,9,by=1),seq(10,30,by=2))
  sft = pickSoftThreshold(data= t(datExpr), networkType = "signed", corFnc="bicor",verbose=2,powerVector=powers)
  save(file=paste("./",acronym_labels[i],"sft.RData",sep = ""),sft)
  
  pdf(file=paste("./",acronym_labels[i],"sft.pdf",sep=""), width=5, height=7)
  par(mfrow=c(2,1))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.9, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
  abline(h=0.8, col="black")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col="black")
  dev.off()
  
  
  ##########################################################################################################
  ##              Run WGCNA when  soft-threshold >=0.8 & sft <= 0.9
  ##########################################################################################################
  
  sft_cutoff <- as.data.frame(sft)
  sft_cutoff_st <- sft_cutoff[which(sft_cutoff$fitIndices.SFT.R.sq >= 0.8),]
  sft_cutoff_end <- sft_cutoff[which(sft_cutoff$fitIndices.SFT.R.sq >= 0.9),]
  
  st <- sft_cutoff_st$fitIndices.Power[1]
  end <- sft_cutoff_end$fitIndices.Power[1]
  
  for (ss in st:end){
    path <- paste(file="./",acronym_labels[i],"/sft",ss,sep="")
    path2 <- paste(file="./",acronym_labels[i],"/sft",ss,"/modSize25",sep="")
    path3 <- paste(file="./",acronym_labels[i],"/sft",ss,"/modSize25/GO_wbackground",sep="")
    path4 <- paste(file="./",acronym_labels[i],"/sft",ss,"/modSize25/Module",sep="")
    
    dir.create(file.path(path))
    dir.create(file.path(path2))
    dir.create(file.path(path3))  
    dir.create(file.path(path4))
    
    
    # calculate TOM
    net = blockwiseModules(datExpr=t(datExpr), maxBlockSize=22000,networkType="signed",corType="bicor", power = ss,
                           mergeCutHeight= 0.1, minModuleSize= 25, pamStage=TRUE, reassignThreshold=1e-6, 
                           saveTOMFileBase=paste(file=path2,"/TOM_signed",sep = ""), saveTOMs=TRUE, verbose = 5, deepSplit=2,nThreads=6)
    save(file=paste(file = path2,"/net_sft",ss,"_modSize25.RData",sep=""),net)
    table(net$colors)
    
    
    # recut the dendrogram manually
    ds = 2; minModSize = 25; dthresh = 0.1; pam = FALSE
    load(paste(file=path2,"/TOM_signed-block.1.RData",sep=""))
    networks=list()
    networks$datExpr=datExpr
    networks$tree = hclust(1-TOM, method="average")
    networks$cut = cutreeHybrid(dendro = networks$tree, pamStage=pam, minClusterSize= minModSize, cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-TOM))
    networks$merged=merged = mergeCloseModules(exprData= t(networks$datExpr), colors = networks$cut$labels, cutHeight=dthresh)
    networks$MEs = moduleEigengenes(t(networks$datExpr), colors=networks$merged$colors, softPower=ss)
    networks$kMEtable = signedKME(t(networks$datExpr), datME = networks$MEs$eigengenes,corFnc = "bicor")
    
    pdf(paste(path2,"/ds2_sft",ss,"_modSize25_dendrogram.pdf",sep=""), width=7.5, height=3.5)
    plotDendroAndColors(networks$tree, colors=labels2colors(networks$merged$colors), dendroLabels = F)
    dev.off()
    
    save(file=paste(file=path2,"/ds2.sft",ss,"_mod25.RData",sep = ""), networks,metadata_tmp) 
    
    ##########################################################################################################
    ##             module-time Point association
    ##########################################################################################################
    datExpr=networks$datExpr
    geneTree=networks$tree
    merged=networks$merged
    modules=merged$colors
    print(length(unique(modules))-1)
    MEs=networks$MEs
    kMEtable=networks$kMEtable
    table(rownames(metadata_tmp) == colnames(datExpr))
    metadata_tmp$age_months=factor(metadata_tmp$age_months,levels = c("3","12","15","18","21","26","28"))
    
    modTrait=data.frame()
    for(i in 2:length(unique(modules))) {
      me = MEs$eigengenes[,i]
      moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
      moduleColor = labels2colors(moduleNumber)
      s = tryCatch(summary(lm(me ~ age_months, data=metadata_tmp))$coefficients,error=function(e){NA})
      if (!is.na(s)){
        for(grp in c("12","15","18","21","26","28")) {
          rowID = paste0("age_months", grp)
          modTrait = rbind(modTrait,
                           data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                                      beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))
        }
      }
    }
    modTrait$fdr=1
    
    modTrait$fdr[modTrait$Group=="12"] = p.adjust(modTrait$p[modTrait$Group=="12"], method = "fdr")
    modTrait$fdr[modTrait$Group=="15"] = p.adjust(modTrait$p[modTrait$Group=="15"], method = "fdr")
    modTrait$fdr[modTrait$Group=="18"] = p.adjust(modTrait$p[modTrait$Group=="18"], method = "fdr")
    modTrait$fdr[modTrait$Group=="21"] = p.adjust(modTrait$p[modTrait$Group=="21"], method = "fdr")
    modTrait$fdr[modTrait$Group=="26"] = p.adjust(modTrait$p[modTrait$Group=="26"], method = "fdr")
    modTrait$fdr[modTrait$Group=="28"] = p.adjust(modTrait$p[modTrait$Group=="28"], method = "fdr")
    
    modTrait$signedLog10fdr = -log10(modTrait$fdr) * sign(modTrait$beta)
    modTrait$signedLog10fdr[modTrait$fdr > .05] = 0
    modTrait$text = signif(modTrait$beta, 1)
    modTrait$text[modTrait$fdr > 0.05] = ""
    modTrait$Module=factor(modTrait$Module,levels = unique(modTrait$Module))
    modTrait$ModuleID <- paste("M",modTrait$moduleNumber,":",modTrait$Module,sep = "")
    
    pdf(file=paste(path2,"/ds2_sft",ss,"_modSize25_SignificantModules.pdf",sep=""), width=13, height=3.5)
    
    g_modules <- ggplot(modTrait, aes(x=ModuleID,y=Group, label=text)) +
      geom_tile(aes(fill=signedLog10fdr),color="grey60") + 
      scale_fill_gradient2(low = "blue", mid = "white", high = "red","[beta]\nsigned\n-log10FDR\n")+
      geom_text(aes(label=text),size=3)+
      theme(axis.text.x = element_text(angle = 325, vjust = 0.2, hjust=0))+
      theme(axis.text = element_text(face="bold",color="black")) + ylab ("Age (months) \n")
    plot(g_modules)
    dev.off()
    
    ##########################################################################################################
    ##            Annotation
    ##########################################################################################################
    ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    Mouse = listAttributes(ensembl)
    filter = "ensembl_gene_id"
    attributes = c("ensembl_gene_id","gene_biotype","external_gene_name","ensembl_transcript_id","external_transcript_name","description")
    annotation = getBM(attributes, mart=ensembl)
    
    # module membership
    member=labels2colors(modules)
    membership=data.frame(id=rownames(datExpr),module=member)
    membership$ensemblID <- gene_ensembl$gene_id[match(membership$id,gene_ensembl$gene_symbol)]
    
    datExpr=networks$datExpr
    genes=rownames(datExpr)
    genes <- as.data.frame(genes)
    genes$EnsemblID <- gene_ensembl$gene_id[match(genes$genes,gene_ensembl$gene_symbol)]
    genes <- as.vector(genes$EnsemblID)
    merged=networks$merged
    modules=merged$colors
    MEs = networks$MEs
    kMEtable=networks$kMEtable
    
    exp_genes <- genes
    for(i in 2:length(unique(modules))){
      
      moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
      moduleColor = labels2colors(moduleNumber)
      moduleGenes = genes[modules==moduleNumber]
      query <- as.vector(moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)])
      
      go = gprofiler2::gost(query,organism = "mmusculus", custom_bg = exp_genes,
                            sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC","TF","MIRNA","CORUM","HP","WP"), correction_method = "fdr",
                            significant = T,ordered_query = T,evcodes=TRUE)
      
      go_results <- as.data.frame(go$result)
      if(nrow(go_results)>0){
        go_results = go_results[order(go_results$p_value),]
        go_results = go_results[go_results$p_value < 0.1,]
        go_results <- apply(go_results,2,as.character)
        write.csv(go_results,file=paste(path2,"/GO_wbackground/go.M",moduleNumber,"_",moduleColor,".csv",sep=""))
      }
    }
    
    ##########################################################################################################
    ##            Cell-type enrichment
    ##########################################################################################################
    pSI.zhang = read.csv(file="Zhang_Human_Neuro2015_pSI_GSE21653.csv",row.names=1,stringsAsFactors = F)
    pSI.zeisel = read.csv(file="Zeisel_level1_Mouse_Science2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)
    pSI.goldmann = read.csv(file="Goldman_levelHybrid_Mouse_NatImmunol2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)
    
    modCellType = data.frame()
    annot <- as.data.frame(rownames(counts_tmp))
    colnames(annot)[1] <- "gene_id"
    annot$EnsemblId <- annotation$ensembl_gene_id[match(annot$gene_id,annotation$external_gene_name)]
    
    modCellType = data.frame()
    
    for(i in 2:length(unique(modules))) {
      print(i)
      me = MEs$eigengenes[,i]
      moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
      moduleColor = labels2colors(moduleNumber)
      genesInMod = as.character(annot$EnsemblId[modules==moduleNumber])
      
      convert_genes <- gorth(genesInMod, source_organism = "mmusculus",
                             target_organism = "hsapiens",  numeric_ns = "",
                             mthreshold = Inf, filter_na = T)
      
      genesInModd <- convert_genes$ortholog_ensg
      
      f.zhang = fisher.iteration(pSI.zhang, genesInModd,p.adjust = F)
      f.zeisel = fisher.iteration(pSI.zeisel, genesInModd,p.adjust = F)
      f.goldman = fisher.iteration(pSI.goldmann, genesInModd,p.adjust = F)
      
      
      modCellType = rbind(modCellType, 
                          data.frame(Dataset="Zhang", Module=moduleColor, moduleNumber= moduleNumber, CellType=rownames(f.zhang), p=f.zhang[,2]),
                          data.frame(Dataset="Zeisel", Module=moduleColor, moduleNumber= moduleNumber, CellType=rownames(f.zeisel), p=f.zeisel[,2]),
                          data.frame(Dataset="Goldman", Module=moduleColor, moduleNumber= moduleNumber, CellType=rownames(f.goldman), p=f.goldman[,2])
      )
      
    }
    
    
    modCellType$fdr = p.adjust(modCellType$p, method = "fdr")
    modCellType$log10fdr = -log10(modCellType$fdr) 
    modCellType$text = signif(modCellType$log10fdr,2)
    modCellType$text[modCellType$fdr > 0.05] = ""
    modCellType$CellType = gsub("Oligodendrocyte", "Oligo", modCellType$CellType)
    modCellType$Module_ID <- paste("M",modCellType$moduleNumber,":", modCellType$Module,sep = "")
    
  }
}


