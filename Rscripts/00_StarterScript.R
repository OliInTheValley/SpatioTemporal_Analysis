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

#A few simple helper functions that will be needed to be loaded before commencing any of the analy
#"Not in" function
'%ni%' <- Negate('%in%')

#Returns all elements of a vector x except the ones in y
everything_but <- function(x,y) {
  return(x[x %ni% y])
}

