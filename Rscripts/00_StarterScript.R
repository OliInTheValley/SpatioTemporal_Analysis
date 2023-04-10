#This script loads required libraries that are used in the following analyses, as well as a few helper functions
#We recommend to run this script before running any of the other scripts to have all the relevant packages at hand


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

#A few simple helper functions
'%ni%' <- Negate('%in%')

everything_but <- function(x,y) {
  return(x[x %ni% y])
}

