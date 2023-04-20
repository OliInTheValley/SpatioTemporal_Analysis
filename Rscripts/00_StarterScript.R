library(readr)
library(DESeq2)
library(tidyr)
library(dplyr)
library(reshape2)
library(gplots)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(VISION)
library(ComplexUpset)
library(lsmeans)
library(ggrepel)
library(scales)
library(Vennerable)
library(rstatix)
library(RColorBrewer)
library(GEOquery)
library(limma)

#A few simple helper functions
'%ni%' <- Negate('%in%')

everything_but <- function(x,y) {
  return(x[x %ni% y])
}





stars.return <- function(input){
  if (is.na(input)){
    return("")
  } else  if (input <= 0.001){
    return("***")
  } else if (input <= 0.01){
    return("**")
  } else if (input <= 0.05){
    return("*")
  } else if (input <= 0.1){
    return(".")
  } else {return("")}
}
