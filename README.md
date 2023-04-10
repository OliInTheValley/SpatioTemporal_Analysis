# Spatio-temporal brain aging map analysis
This GitHub repository documents the core analyses of a spatiotemporal RNA-seq study on mouse brain aging, including 1,076 samples from 15 regions spanning 7 ages and two rejuvenating interventions. 
The analyses were performed within an R 3.6 enviornment.

We recommend to re-create the directory setup provided here, consisting of:   
1. Rscripts - a folder where the core scripts are stored
2. input_data - a folder containing raw counttables, metadata etc. We provide here only the bulkseq data, as single-cell/spatial data exceed the size limit
3. R_objects - a folder where processed R objects, such as Deseq2 objects are stored so they can be loaded again by later analyses
4. Signature_repository - a folder where we store VISION score objects in the form of R.bin files that can be loaded when needed
5. Output_tables - at some points in the scripts, we create e.g. lists of genes. These sould be stored here  


# Sofware enviornment
R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X  13.2

The following packages and versions are required and should be loaded prior to running the analysis    
Biobase_2.46.0  
BiocGenerics_0.32.0  
BiocParallel_1.20.1  
ComplexUpset_1.3.3  
DelayedArray_0.12.3  
DESeq2_1.26.0  
dplyr_1.0.9  
emmeans_1.7.5  
GenomeInfoDb_1.22.1  
GenomicRanges_1.38.0  
ggplot2_3.3.6  
ggrepel_0.9.1  
gplots_3.1.3  
IRanges_2.20.2  
lsmeans_2.30-0  
matrixStats_0.62.0  
pheatmap_1.0.12  
readr_2.1.2  
reshape2_1.4.4  
rstatix_0.7.0  
S4Vectors_0.24.4  
scales_1.2.0  
Seurat_3.1.2  
SummarizedExperiment_1.16.1  
tidyr_1.2.0  
Vennerable_3.1.0.9000  
VISION_3.0.0  
