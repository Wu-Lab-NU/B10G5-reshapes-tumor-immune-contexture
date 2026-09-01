#load libraries and set directories 
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

setwd("/Users/jem/Documents/Manuscript_Source_Codes") #update path as needed 

NK_data <- readRDS("NK_data.rds") #the original NK object 
NK_data <- NK_data[!grepl("^Rp[sl]", rownames(NK_data)) & !grepl("^mt-", rownames(NK_data)),] # remove mitochondrial and ribosomal genes 
head(NK_data@meta.data) #inspect - should contain "NK.subset" and "Tx" in metadata
p0 <- FeaturePlot(NK_data, features = c('Cd3d','Cd3e','Trac','Kdm6b','Tcf7','Itgae'), reduction = "umap") #note we can see a small CD3d/e+ cluster

#re-do batch correction for NK
NK.batch <- NK_data

NK.batch <- FindVariableFeatures(NK.batch, selection.method = "vst", nfeatures = 2000)
NK.batch <- ScaleData(NK.batch, features = rownames(NK.batch), verbose = FALSE)
NK.batch <- RunPCA(NK.batch, npcs = 50, verbose = F) #Warning: Number of dimensions changing from 50 to 30 - changed to 50
NK.batch <- RunHarmony(NK.batch, "orig.ident", assay.use="RNA", plot_convergence = TRUE, max.iter.harmony = 20)
NK.batch <- RunUMAP(NK.batch, dims = 1:15, reduction = "harmony") #originally set at 30
NK.batch <- FindNeighbors(NK.batch, reduction = "harmony", dims = 1:15) #originally set at 30
NK.batch <- FindClusters(NK.batch, algorithm =2, resolution = 0.3) #originally set at 0.4

p2 <- DimPlot(
  NK.batch,
  reduction = "umap",
  group.by = c("orig.ident", "Tx", "seurat_clusters"),
  combine = T, label.size = 2
)


p3 <- DotPlot(NK.batch, features = c('Kdm6b',"Tcf7","Sell", "Klf2", "S1pr1", "Zeb2", "Eomes", "Gzma", "Prf1", 
                                     "Il18rap","Il18r1", "Irf8", "Klra4", 'Itga1','Ntan1', 'Ltb','Cxcr6','Gpr34', 'Itgae', 'Gzmc', 'Cd200r4','Cd3d','Cd3e','Trac')) + RotatedAxis()

p4 <- FeaturePlot(NK.batch, features = c('Cd3d','Cd3e','Trac','Kdm6b','Tcf7','Itgae'), reduction = "umap")

(p2 / p3 )| p4 #cluster 3 appears to be T cells as it is enriched in CD3d/e and Trac

NK_fate <- subset(NK.batch, idents = c('3'), invert = T) #remove cluster 3
head(NK_fate@meta.data)
levels(NK_fate$seurat_clusters) #"0" "1" "2"
p5 <- FeaturePlot(NK_fate, features = c('Cd3d','Cd3e','Trac','Kdm6b','Tcf7','Itgae'), reduction = "umap")

#update NK.subset
NK_fate@meta.data <- NK_fate@meta.data %>% #Assign new ident "NK.subset"
  mutate(NK.subset = case_when(
    seurat_clusters == "0" ~ "trNK",
    seurat_clusters == "1" ~ "cytoNK",
    seurat_clusters == "2" ~ "stemNK",
    TRUE ~ seurat_clusters # If there are other values in orig.ident, keep them unchanged
  ))
head(NK_fate@meta.data) # Verify the new metadata column

Idents(NK_fate) <- 'NK.subset'
my_level = c('trNK',"stemNK","cytoNK")
NK_fate$NK.subset <- factor(NK_fate$NK.subset, levels = my_level) #using factor() solved the problem the level() cannot
levels(NK_fate) #check

#save rds object for batch corrected
saveRDS(NK_fate, file = "NK_data_batch.rds")
