#load libraries and set directories ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(harmony)
library(slingshot)
library(viridis)
library(openxlsx)
library(presto) #devtools::install_github('immunogenomics/presto')
library(readr) #install.packages("readr")
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")  # Just in case, prevent OpenMP crash
library(SCP) 

setwd("/Users/jem/Documents/Manuscript_Source_Codes") #update path as needed 
NK_fate <- readRDS("NK_data_batch.rds") #read NK object

#Trajectory Inference with Slingshot (Figure 3D-E & Supp Figure 2)----
##Figure 3D
test_1 <- RunSlingshot(srt = NK_fate, group.by = "NK.subset", reduction = "umap", 
                       start = c('stemNK'), 
                       end = c('trNK','cytoNK'), 
                       align_start = T, extend = "n", stretch = 0) #extend = "n" originally
ggsave(filename = "Output/Re-Batch/600dpi-final/FigS2-NK-slingshot.png", device = 'png', height = 5, width = 5, dpi = 600) 

##Figure 3E & FigS2
DynamicPlot(
  srt = test_1, lineages = c("Lineage1", "Lineage2"), group.by = "NK.subset", 
  features = c("Cxcr6", 'Itga1',"Itgae","Prf1",'Gzma','Tcf7','Kdm6b'),
  compare_lineages = TRUE, compare_features = FALSE
)
ggsave(filename = "Output/Re-Batch/600dpi-final/FigS2-NK-dynamic-plot.png", device = 'png', height = 6, width = 9, dpi = 600) #Figure S2 + Figure 3E

DynamicPlot(
  srt = test_1, lineages = c("Lineage1", "Lineage2"), group.by = "NK.subset", 
  features = c('Tcf7','Kdm6b'),
  compare_lineages = TRUE, compare_features = FALSE, 
  line.size = 0.5, pt.size = 0.5, 
)
ggsave(filename = "Output/Re-Batch/600dpi-final/Fig3e-NK-dynamic-plot_Tcf1+Kdm6b.png", device = 'png', height = 3, width = 9, dpi = 600) #Figure S2 + Figure 3E

DynamicPlot(
  srt = test_1, lineages = c("Lineage1", "Lineage2"), group.by = "NK.subset", 
  features = c("Cxcr6", 'Itga1',"Itgae","Prf1",'Gzma','Tcf7'),
  compare_lineages = TRUE, compare_features = FALSE
)
ggsave(filename = "Output/Re-Batch/600dpi-final/Fig3e+S2-NK-dynamic-plot-size2.png", device = 'png', height = 5, width = 10, dpi = 600) #Figure S2 + Figure 3E
