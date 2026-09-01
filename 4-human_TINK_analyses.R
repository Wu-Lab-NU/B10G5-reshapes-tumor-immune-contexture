#R codes for figure 3 - analyses human pan-cancer NK scRNA-seq dataset 
#For figure 3H-J + Figure S6C-D
#Compiled and last updated on - 2025-07-24

# Load the required libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(openxlsx)

#Construct Seurat Object ----
#raw data (GSE212890) downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212890
setwd("/Users/jem/Documents/Manuscript_Source_Codes/human_TINK") #update path as needed 
#setwd("/Users/jem/Documents/Manuscript_Source_Codes") 

# Define the paths to the files
barcodes_file <- "raw/NK_rawcount.barcodes"
genes_file <- "raw/NK_rawcount.genes"
counts_file <- "raw/NK_rawcount.mtx"
metadata_file <- "raw/GSE212890_NK_metadata.csv"

# Read in the data
barcodes <- read.csv(gzfile(barcodes_file), header = FALSE)
genes <- read.csv(gzfile(genes_file), header = FALSE)
counts <- readMM(gzfile(counts_file))
metadata <- read.csv(metadata_file,row.names = 1)
head(counts)

# Check the dimensions
print(dim(counts))  # Should match the number of genes and barcodes
print(length(barcodes$V2))
print(length(genes$V2),row.names = 1)

#remove row 1 from both barcodes and genes which are null
barcodes <- barcodes[-1,]
genes <-genes[-1,]

# Rename the columns and rows appropriately
colnames(counts) <- genes$V2
rownames(counts) <- barcodes$V2
transposed_counts <- t(counts)

# Create a Seurat object
TINK <- CreateSeuratObject(counts = transposed_counts)

# Add metadata to the Seurat object
TINK <- AddMetaData(TINK, metadata)

# Save the Seurat object (optional)
saveRDS(TINK, file = "TINK-unprocessed.rds")


#QC + clustering----

#Standard pre-processing workflow with Seurat
TINK[["percent.mt"]] <- PercentageFeatureSet(TINK, pattern = "^MT-")
TINK[['percent.ribo']] <- PercentageFeatureSet(TINK, pattern = "^RP[SL]")

head(TINK@meta.data, 10)
#TINK <- readRDS("TINK-unprocessed.rds")


# Visualize QC metrics as a violin plot
VlnPlot(TINK, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 4) 
plot1 <- FeatureScatter(TINK, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TINK, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
TINK <- subset(TINK, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #check different percent.mt cut-offs


#Normalizing the data
TINK <- NormalizeData(TINK)

TINK <- FindVariableFeatures(TINK, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(TINK), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(TINK)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(TINK)
TINK <- ScaleData(TINK, features = all.genes)

#Perform linear dimensional reduction
TINK <- RunPCA(TINK, features = VariableFeatures(object = TINK))
VizDimLoadings(TINK, dims = 1:2, reduction = "pca")

ElbowPlot(TINK,ndims = 50) 

#Cluster the cells 
TINK <- FindNeighbors(TINK, dims = 1:10) #use dim = 10 default setting 
TINK <- FindClusters(TINK, resolution = 0.32) #0.32 looks the best - 0.35 = 0.4 no diff

head(Idents(TINK), 5) 
TINK <- RunUMAP(TINK, dims = 1:10)
DimPlot(TINK, reduction = "umap")
DimPlot(TINK, reduction = "umap",split.by = 'meta_tissue')

VlnPlot(TINK, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 4) 
saveRDS(TINK, file = "TINK.rds")
#TINK <- readRDS("TINK.rds")


#check QC criteria post clustering + final clean-up ----
VlnPlot(TINK, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 4) 
table(TINK$seurat_clusters)
#cluster 5 appear to have low nFeature_RNA and nCount_RNA
#cluster 6 has < 100 cells
plot1 <- FeatureScatter(TINK, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TINK, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

VlnPlot(TINK, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 4) 

TINK <- subset(TINK, idents = c('5','6'), invert = TRUE) #remove cluster 5 & 6  - low quality AND/OR low cell number
levels(TINK)
DimPlot(TINK, reduction = "umap")
DimPlot(TINK, reduction = "umap",split.by = 'meta_tissue')

#subset tumor-infiltrating NK cells ----
tumor <- subset(x = TINK, subset = meta_tissue == 'Tumor') #meta_tissue contains "tumor" and matching "normal" tissue

#Cluster Annotation ----
p1 <- DimPlot(tumor, reduction = "umap")
ggsave(p1, filename = "output/0-umap_seurat_clusters_tumor.pdf", device = 'pdf', height = 3, width = 3)

markers <- FindAllMarkers(tumor, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
ph <- DoHeatmap(tumor, features = top15$gene) + NoLegend()
ggsave(ph, filename = "output/0-DoHeatmap_top15_log2Fc.png", device = 'png', height = 12, width = 8.5)

markers <- markers[markers$p_val_adj <= 0.05,]
write.xlsx(markers, file = 'output/0-markers_cluster_tumorOnly.xlsx',rowNames = T)

#annotate cells based on key markers
markers <- c('NCAM1','FCGR3A','NCR1','CXCR6','ITGA1','ITGAE','TCF7','EOMES','TBX21','PRF1','GZMA','GZMB','MKI67')
ph <- VlnPlot(tumor, features = markers,ncol = 3)
ggsave(ph, filename = "output/0-NK_subsetMarkers_Violin.pdf", device = 'pdf', height = 11, width = 8.5)

ph <- DotPlot(tumor, features = markers, dot.scale = 6,cols = c("blue", "red")) + RotatedAxis() + coord_flip()
ggsave(ph, filename = "output/0-DotPlot_NK_Subseet_Markers.pdf", device = 'pdf', height = 4, width = 5)

levels(tumor$seurat_clusters)
tumor@meta.data <- tumor@meta.data %>%
  mutate(NK.subset = case_when(
    seurat_clusters == "0" ~ "cytoNK", #CD16+Perforin-hi
    seurat_clusters == "1" ~ "cytoNK",
    seurat_clusters == "2" ~ "cytoNK",
    seurat_clusters == "3" ~ "trNK", #CXCR6+CD49a+CD103+
    seurat_clusters == "4" ~ "stemNK", #TCF7-hi
    TRUE ~ seurat_clusters # If there are other values in orig.ident, keep them unchanged
  ))
head(tumor@meta.data) # Verify the new metadata column
saveRDS(tumor, file = "TINK_annotated.rds") 

#tumor <- readRDS(file = "/Users/jem/Documents/Manuscript_Source_Codes/human_TINK/TINK_annotated.rds") 
Idents(tumor) <- "NK.subset"
levels(tumor) 

plot <- DimPlot(tumor, reduction = "umap", label = F, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 4))) #4x5
ggsave(plot, filename = "output/Fig3/Fig3h-umap-NK-subset.pdf", device = 'pdf', height = 4, width = 5) #(Figure 3H)
ggsave(plot, filename = "output/Fig3/Fig3h-umap-NK-subset.png", device = 'png', height = 4, width = 5, dpi = 300) #(Figure 3H)

#Redo FindMarkers for the annotated clusters
markers <- FindAllMarkers(tumor, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15
ph <- DoHeatmap(tumor, features = top15$gene) + NoLegend()
ggsave(ph, filename = "output/Fig3/FigS6c-Annotated_DoHeatmap_top15_log2Fc.png", device = 'png', height = 12, width = 8.5, dpi = 300)  #(Matching Supp Figure A)
ggsave(ph, filename = "output/Fig3/FigS6c-Annotated_DoHeatmap_top15_log2Fc.pdf", device = 'pdf', height = 12, width = 8.5)  #(Matching Supp Figure A)

markers <- markers[markers$p_val_adj <= 0.05,]
write.xlsx(markers, file = 'output/1-markers_cluster_tumorOnly+annotated.xlsx',rowNames = T)


#Select Markers
markers <- c('NCAM1','FCGR3A','NCR1','CXCR6','ITGA1','ITGAE','TCF7','EOMES','TBX21','PRF1','GZMA','GZMB','MKI67')
ph <- VlnPlot(tumor, features = markers,ncol = 3)
ggsave(ph, filename = "output/1-NK_subsetMarkers_Violin-merged.pdf", device = 'pdf', height = 11, width = 8.5)

ph <- DotPlot(tumor, features = markers, dot.scale = 6,cols = c("blue", "red")) + RotatedAxis() + coord_flip()
ggsave(ph, filename = "output/1-DotPlot_NK_Subset_Markers-merged.pdf", device = 'pdf', height = 4, width = 5)


library(scales)
identities <- levels(tumor@active.ident)
my_cols <- hue_pal()(length(identities))
show_col(my_cols)
library(scCustomize)
main <- c('TCF7','PRF1','FCGR3A','ITGAE','CXCR6')
ph <- Stacked_VlnPlot(seurat_object = tumor, features = main, x_lab_rotate = TRUE, colors_use = my_cols, pt.size = 0.01, alpha = 0.4)
ggsave(ph, filename = "output/Fig3/Fig3i-Stacked-Violin-NK-subset-marker.pdf", device = 'pdf', height = 5, width = 3.5) #(Figure 3I)
ggsave(ph, filename = "output/Fig3/Fig3i-Stacked-Violin-NK-subset-marker.png", device = 'png', height = 5, width = 3.5, dpi = 300) #(Figure 3I)

#dotplot for stemNK&B10G5-up overlap (identified with mouse scRNA-seq) ----
genes <- toupper(c('Npm1','Satb1','Chd7','Chd4','Chd3','Baz1a','Emsy','Smarca4','Nasp','Myb','Bptf','Kdm6b','Kdm2b','Ash1l'))
ph <- DotPlot(tumor, features = rev(genes), dot.scale = 6,cols = c("blue", "red")) + RotatedAxis() + coord_flip() 
ggsave(ph, filename = "output/1-DotPlot_stemNK&B10G5-overlap-human-TINK.pdf", device = 'pdf', height = 4, width = 5) #(Figure 3J)

#dotplot for SATB1. KDM6B and CHD4 only 
genes <- c("SATB1","KDM6B","CHD4")
ph <- DotPlot(tumor, features = rev(genes), dot.scale = 6,cols = c("blue", "red")) + RotatedAxis() + coord_flip() 
ggsave(ph, filename = "output/1-DotPlot_stemNK&B10G5-SATB1+KDM6B+CHD4.pdf", device = 'pdf', height = 4, width = 5) #(Figure 3J - alternative)


#Exploration only ----
genes <- c('LAG3','TOX','XCL1','CCL5','FASLG','TNFSF10','MYC','LDHA','PGAM1')
genes <- c('SLC7A8','ALDOB','PGM5','SLC2A2')
genes <- c('TOX','DGAT1','FASN') 
genes <- c('LAG3','TOX','TCF7', 'KLRE1','KLRG1','KLRD1','KLRK1')

#beta-catenin
genes <- c('CTNNB1','TCF7') 
ph <- DotPlot(tumor, features = genes, dot.scale = 6,cols = c("blue", "red")) + RotatedAxis() + coord_flip() 
p <- VlnPlot(tumor, features = c('CTNNB1'))


#for thesis-KDM6B-hi vs. lo ----
ph <- VlnPlot(tumor, features = 'KDM6B')
ggsave(ph, filename = "Output/thesis/Violin-Kdm6b-tumor.pdf", device = 'pdf', height = 3, width = 4)

Kdm6b_high <- WhichCells(object = tumor, expression = KDM6B >= 1.5) 
Kdm6b_low <- WhichCells(object = tumor, expression = KDM6B < 1.5) 
Kdm6b_marker <- FindMarkers(tumor, ident.1 = Kdm6b_high, ident.2 = Kdm6b_low, verbose = FALSE)
Kdm6b_marker <- Kdm6b_marker[Kdm6b_marker$p_val_adj <= 0.05,] 

#save the DEG 
write.csv(Kdm6b_marker,file= "Output/thesis/Kdm6b-hi_vs_lo_markers.csv",quote=F)

#volcano plot-----
tfnames_show <- c('REL','NFKB1','RELB','NFKB2','NFKBIZ','NFKBIA','NFKBID')
Kdm6b_marker <- Kdm6b_marker %>%
  mutate(text_tf = ifelse(rownames(Kdm6b_marker) %in% tfnames_show, rownames(Kdm6b_marker), NA))
Kdm6b_marker <- Kdm6b_marker[-1,] #remove Kdm6b from the list

#add dashed reference line for y = -log10(p_val_adj) = 1, x = avg_log2FC = 1; points below should be grey
Kdm6b_marker$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
Kdm6b_marker$diffexpressed[Kdm6b_marker$avg_log2FC > 1 & Kdm6b_marker$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.05, set as "DOWN"
Kdm6b_marker$diffexpressed[Kdm6b_marker$avg_log2FC < -1 & Kdm6b_marker$p_val_adj < 0.05] <- "DOWN"

library(ggrepel)

p1 = ggplot(Kdm6b_marker, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed)) +
  labs(x = "avg_log2FC", y = "-Log10(P.adj)") +
  theme_classic() +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(alpha = 0.8, shape = 16) + 
  geom_text_repel(aes(label=text_tf), size = 4, max.overlaps = 20) +
  scale_x_continuous(breaks = c(seq(-4, 4, 1)), # Modify x-axis tick intervals  
                     limits = c(-4, 4))+
  theme(
    axis.title.x = element_text(size = 12),  # Set x-axis label font size and style
    axis.title.y = element_text(size = 12)   # Set y-axis label font size and style
  ) +
  scale_color_manual(values = c("red", "grey", "blue"), 
                     labels = c("Kdm6b-lo", " ", "Kdm6b-hi")) 
ggsave(p1,filename = "output/volcano_Kdm6b-hi_v_lo.pdf", device = 'pdf', height = 4, width = 5 )


