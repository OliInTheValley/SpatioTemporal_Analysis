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

setwd("./WGCNA/")

library(gprofiler2)
library(WGCNA)
library(ggplot2)
library(readr)
library(DESeq2)
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

