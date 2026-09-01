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

#Figure 1 and related supps ----
my_col = c('#F68282','#619CFF',"#0CB702")

Idents(NK_fate) <- 'seurat_clusters'
levels(NK_fate) #check

ph0 <- DimPlot(
  NK_fate,
  reduction = "umap",
  group.by = c("Tx", "seurat_clusters","NK.subset"),
  combine = T, label.size = 2
) #for preview only

##Figure 1B
Idents(NK_fate) <- 'seurat_clusters'

ph1 <- DimPlot(NK_fate, reduction = "umap", label = T, label.size = 4.5, cols = my_col) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "none", axis.title = element_text(size = 18), legend.text = element_text(size = 18),axis.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) #4x5
ggsave(ph1, filename = "Output/Re-Batch/600dpi-final/Fig1b-NK_UMAP-all.png", device = 'png', height = 4, width = 4.5, dpi = 600) #Figure 1B

ph2 <- DimPlot(NK_fate, reduction = "umap", label = TRUE, label.size = 0, group.by = c("Tx"), cols = c("#7a7a7a","#cb5c5c"), alpha = 0.3) + xlab("UMAP 1") + ylab("UMAP 2")  +
    theme(legend.position = "none", plot.title = element_blank(), axis.title = element_text(size = 18), legend.text = element_text(size = 18),axis.text = element_text(size = 18)) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) #4x5

ggsave(ph2, filename = "Output/Re-Batch/600dpi-final/Fig1b-NK_UMAP-overlaid.png", device = 'png', height = 4, width = 4.5, dpi = 600) #Figure 1B - overlaid

##Figure 1D-E 
ph3 <- DimPlot(NK_fate, reduction = "umap", label = T, label.size = 4.5, cols = my_col, split.by = 'Tx') + xlab("UMAP 1") + ylab("UMAP 2") + 
  guides(colour = guide_legend(override.aes = list(size = 4))) #4x7
ggsave(ph3, filename = "Output/Re-Batch/600dpi-final/Fig1d-e_dimplot_untreated+B10G5.png", device = 'png', height = 4, width = 7, dpi = 600) #Figure 1D-E 

ph3a <- DimPlot(
  subset(NK_fate, subset = Tx == "No Tx"),
  reduction = "umap",
  label = TRUE,
  label.size = 4.5,
  cols = my_col
) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("No Tx") +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  NoLegend()

ph3b <- DimPlot(
  subset(NK_fate, subset = Tx == "B10G5"),
  reduction = "umap",
  label = TRUE,
  label.size = 4.5,
  cols = my_col
) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  ggtitle("B10G5") +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  NoLegend()

ph3 <- ph3a + ph3b
ggsave(ph3, filename = "Output/Re-Batch/600dpi-final/Fig1d-e_dimplot_untreated+B10G5.png", device = 'png', height = 4, width = 7, dpi = 600) #Figure 1D-E 

##Figure 1C
select.top20 <- c('Cxcr6','Itga1','Itgae','Rgs1','Xcl1','Tcf7','Eomes','Sell','Myb','Satb1','Prf1','Gzma','Irf8','Klf2','Klrg1')
ph4 <- DotPlot(NK_fate, features = rev(select.top20), dot.scale = 6) + RotatedAxis() + coord_flip() 
ggsave(ph4, filename = "Output/Re-Batch/600dpi-final/Fig1c-DotPlots_core_markers_revision.png", device = 'png', height = 4, width = 4.5, dpi = 600) #Figure 1C

##Figure 1F
count <- table(Idents(NK_fate),NK_fate$Tx) #Idents(NK_fate) <- 'NK.subset'
write.xlsx(count, file = "Output/Re-Batch/cell_count.xlsx", rowNames = T) #For Figure 1F

##Figure 1G
Idents(NK_fate) <- 'NK.subset'
markers <- c('Cxcr6','Itga1','Tcf7', 'Eomes','Prf1')
p1 <- FeatureStatPlot(NK_fate, stat.by = markers, group.by = "NK.subset", palcolor = c('#F68282',"#0CB702",'#619CFF'),add_point = T, stack = T, ylab = '') #graphing function in SCP; base = ggplpt
ggsave(p1, filename = "Output/Re-Batch/600dpi-final/Fig1g-ViolinPlots_subset_markers_validation.png", device = 'png', height = 4, width = 4, dpi = 600) 

##Figure S1C
Idents(NK_fate) <- 'seurat_clusters'
NK.Func.Het <- c("Klrk1","Ncr1","Il2rb","Il18r1","Il18rap","Cd247","Fcer1g","Gzma","Gzmb","Prf1",'Fasl','Tnfsf10',"Ifng","Ccl3","Ccl4",'Ccl5','Xcl1') 
p3 <- VlnPlot(NK_fate, features = NK.Func.Het,ncol = 3,cols = my_col)
ggsave(p3, filename = "Output/Re-Batch/600dpi-final/FigS1C-NK_markers1_Violin_withDeathLigand.png", device = 'png', height = 12, width = 8.5, dpi = 600) #Figure S1C - related to figure 1

##Figure S1D
NK.markers <- c('Cxcr6','Itga1', 'Itgae','Rgs1','Prf1','Gzma','Irf8','Klrg1','Tcf7','Eomes','Myb','Satb1')
p4 <- FeaturePlot(NK_fate, features = NK.markers, ncol = 4) 
ggsave(p4, filename = "Output/Re-Batch/600dpi-final/FigS1D-NK_markers_featureplot.png", device = 'png', height = 8, width = 12, dpi = 600) #Figure S1D - related to figure 1

##Table S1 
#Top20-NK cluster markers 
Idents(NK_fate) <- 'seurat_clusters'
NK.markers <- FindAllMarkers(NK_fate, only.pos = TRUE)
NK.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write.xlsx(top20, file = "Output/Re-Batch/600dpi-final/top20_byClusters.xlsx", rowNames = F) #table S1

#Figure 2 and related supps ----
##Figure S5A
Inhibitory <- c("Lag3","Klrc1","Klrg1",'Havcr2','Pdcd1','Tigit')
p5 <- VlnPlot(NK_fate, features = Inhibitory,cols = my_col, ncol = 3)
ggsave(p5, filename = "Output/Re-Batch/600dpi-final/FigS5A-NK_inhibitory_receptors.png", device = 'png', height = 6, width = 9, dpi = 600) #Figure S5A - related to figure 2

#check percent PD-1+ trNK

##Figure S5B
p6 <- DotPlot(NK_fate, features = Inhibitory,cols = c("blue", "red")) + RotatedAxis() + coord_flip() 
ggsave(p6, filename = "Output/Re-Batch/600dpi-final/FigS5B_dotplot_inhibitory.png", device = 'png', height = 3, width = 5, dpi = 600) #Figure S5B

##Figure 2A-B (violin)
Idents(NK_fate) <- 'Tx'
my_level = c("No Tx","B10G5")
NK_fate$Tx <- factor(NK_fate$Tx, levels = my_level) 
levels(NK_fate) #check

Idents(NK_fate) <- 'NK.subset'
ph <- VlnPlot(NK_fate, features = c("Lag3","Klrc1"),idents = 'trNK', group.by = 'Tx', cols = c("#7a7a7a","#cb5c5c"))
ggsave(ph, filename = "Output/Re-Batch/600dpi-final/Fig2a-b-Violin_trNK-only-Lag3+Krlc1_byTx.png", device = 'png', height = 3, width = 6, dpi = 600) #Figure2A-B

#use SCP to perform comparison and obtain significance value
test <- c("Lag3","Klrc1")
p1 <- FeatureStatPlot(NK_fate, stat.by = test, group.by = "NK.subset", split.by = "Tx", palcolor = c("#7a7a7a","#cb5c5c"), add_point = T, stack = T
                      ,comparisons = TRUE,sig_label = c("p.format")
)
ggsave(p1, filename = "Output/Re-Batch/600dpi-final/Fig2a-b-Violin_Lag3+Klrc1_subset+Tx+p-value.png", device = 'png', height = 6, width = 6) #p-val used for Figure2A-B


##FigS5c - Klrg1 expression pattern 
ph <- VlnPlot(NK_fate, features = c('Klrg1'),idents = 'cytoNK', group.by = 'Tx', cols = c("#7a7a7a","#cb5c5c"))
p1 <- FeatureStatPlot(NK_fate, stat.by = c('Klrg1'), group.by = "NK.subset", split.by = "Tx", palcolor = c("#7a7a7a","#cb5c5c"), add_point = T, stack = T
                      ,comparisons = TRUE,sig_label = c("p.format"))
ggsave(ph, filename = "Output/Re-Batch/600dpi-final/FigS5D-Violin_klrg1_cytoNK_only.png", device = 'png', height = 3, width = 3, dpi = 600) 
ggsave(p1, filename = "Output/Re-Batch/600dpi-final/FigS5D-Violin_klrg1_subset+Tx_sig.png", device = 'png', height = 3, width = 5) 


##Fig2E
#B10G5 vs. no Tx bulk DEG
Idents(NK_fate) <- 'Tx'
markers1 <- FindMarkers(NK_fate, ident.1 = "B10G5", ident.2 = "No Tx", verbose = FALSE)
markers1 <- markers1[markers1$p_val_adj <= 0.05,] #need to filter on p.adj <= 0.05
write.xlsx(markers1,file= "Output/Re-Batch/600dpi-final/B10G5vsNoTx_markers.xlsx",rowNames = T) 
#refer to "3-mouse_TINK_DE_analyses(Tx).R" for pathway analysis

##Fig2F
#Module Scores for oxidative stress/damage signatures
integrated_stress_response <- c('Gclm','Glrx','Gpx1','Prdx5','Prdx6','Phgdh','Taldo1', #anti-oxidant - Nrf2 targets
                                'Psma1','Psma4','Psma7','Psmd6', #proteosome unit - Nrf2 targets
                                #'Hspa5','Calr','Ppdia3','Ppdia6', #ER protein folding - involved but not necessarily transcriptionally regulated
                                'Ddit3','Atf4','Txnip','Osgin1','Gadd45g', 'Dapk2') #stress-induced signaling and apoptosis
NK_fate <-AddModuleScore(NK_fate,features = list(integrated_stress_response), name = 'stress')

#all clusters
ph <- FeatureStatPlot(NK_fate, stat.by = c('stress1'), group.by = "NK.subset",split.by = "Tx", palcolor = c("#7a7a7a","#cb5c5c"),add_point = T, stack = F, ylab = '',
                      ,comparisons = TRUE,sig_label = c("p.format")) 

#just trNK + cytoNK
Idents(NK_fate) <- 'NK.subset'
NK_subset2 <- subset(NK_fate, idents = c('trNK', 'cytoNK'))
NK_subset2$Tx <- factor(NK_subset2$Tx,levels = c("No Tx","B10G5")) 

ph2 <- FeatureStatPlot(NK_subset2, stat.by = c('stress1'), group.by = 'NK.subset', split.by = "Tx",palcolor = c("#7a7a7a","#cb5c5c"), add_point = T, stack = T
                       ,comparisons = TRUE,sig_label = c("p.format"), ylab = '', legend.position = 'none')

ggsave(ph2, filename = "Output/Re-Batch/600dpi-final/Fig2f-integrated_stress_response_subset+Tx.png", device = 'png', height = 2.5, width = 3.2, dpi = 600) #Figure 2f #6

##Figure S6A
ht8 <- GroupHeatmap(NK_subset2,
                    features = integrated_stress_response,
                    group.by = "NK.subset", 
                    split.by = "Tx",
                    cluster_rows = F, cluster_columns = F, cluster_row_slices = F, cluster_column_slices = F,
                    add_dot = TRUE, add_reticle = F, #heatmap_palette = "viridis",
                    nlabel = 0, show_row_names = TRUE,
                    ht_params = list(row_gap = unit(0, "mm"), row_names_gp = gpar(fontsize = 12)),
                    cell_split_palcolor = c("#7a7a7a","#cb5c5c","#7a7a7a","#cb5c5c"), 
                    group_palcolor = list(c('#F68282','#619CFF')),
                    label_size = 8)
png("Output/Re-Batch/600dpi-final/FigS6A-dotplot_integrated_stress_response_trNK+cytoNK-only.png", height = 5, width = 5, bg = 'transparent', res = 600, unit = 'in') 
ht8$plot
dev.off()

#Figure 3A-C and related supps ----
#stemNK subset markers
Idents(NK_fate) <- 'NK.subset' 
NK.markers <- FindAllMarkers(NK_fate, only.pos = TRUE)
NK.markers <- NK.markers[NK.markers$p_val_adj <= 0.05, ]
write.xlsx(NK.markers, file = "Output/Re-Batch/600dpi-final/NK_cluster_markers.xlsx", rowNames = T)

head(NK.markers)
stemNK.markers <- NK.markers[NK.markers$cluster == 'stemNK', ]
stemNK.markers <- stemNK.markers$gene #a vector of characters 
length(stemNK.markers) #196

#obtain B10G5 upregulated DEGs 
Idents(NK_fate) <- 'Tx'
levels(NK_fate)
markers1 <- FindMarkers(NK_fate, ident.1 = "B10G5", ident.2 = "No Tx", verbose = FALSE)
markers1 <- markers1[markers1$p_val_adj <= 0.05,] #already saved
B10G5.up <- markers1[markers1$avg_log2FC > 0,]
B10G5.up <- row.names(B10G5.up)
length(B10G5.up) #565

##Figure 3A
library(ggVennDiagram) #install.packages("ggVennDiagram")
x <- list(A = stemNK.markers,
          B = B10G5.up)

plot1 <- ggVennDiagram(x,category.names = c("stemNK","B10G5 up"), label = 'count') + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + #set custom color
  theme(legend.position = 'none') + 
  coord_flip()
ggsave(plot1, filename = "Output/Re-Batch/600dpi-final/Fig3a-Venn-stemNK&B10G5-up.png", device = 'png', height = 3, width = 4, dpi = 600, bg = 'transparent') #Figure 3A - horizontal 

##Figure 3B
#perform EnrichR pathway analysis for the overlapping genes
overlap <- intersect(stemNK.markers,B10G5.up)
length(overlap) #175
cat(
  "overlap <- c(",
  paste0("", overlap, "", collapse = ","),
  ")\n"
)

library(enrichR)
listEnrichrSites() #connect to online database

#select your favourite databases
dbs <- c("KEGG_2019_Mouse","Reactome_2022","MSigDB_Hallmark_2020","ENCODE_TF_ChIP-seq_2015","GO_Biological_Process_2015")

#perform enrichR 
overlap.enrich = enrichr(overlap, dbs) #a list for the databases selected
#head(overlap.enrich[[1]])
for (j in 1:length(dbs)) {
  term = overlap.enrich[[j]]
  term = subset(term, select = -c(Old.P.value,Old.Adjusted.P.value)) #remove empty columns 
  term = subset(term, Adjusted.P.value <= 0.05) #filter terms with p.adj <= 0.05
  filename = paste("Output/Re-Batch/03-StemNK&B10G5-up-overlap-EnrichR-",dbs[j],".xlsx", sep="")
  write.xlsx(term,filename,rowNames = FALSE)
}
#note: use graphpad to make custom bargraphs...

##Fig3C
genes <- c('Npm1','Satb1','Chd7','Chd4','Chd3','Emsy','Myb','Kdm6b','Kdm2b','Ash1l') #genes identified from Figure 3A-B

ht8 <- GroupHeatmap(NK_fate,
                    features = genes,
                    group.by = "NK.subset", 
                    split.by = "Tx",
                    cluster_rows = F, cluster_columns = F, cluster_row_slices = F, cluster_column_slices = F,
                    add_dot = TRUE, add_reticle = F, #heatmap_palette = "viridis",
                    nlabel = 0, show_row_names = TRUE,
                    ht_params = list(row_gap = unit(0, "mm"), row_names_gp = gpar(fontsize = 12)),
                    cell_split_palcolor = c("#7a7a7a","#cb5c5c","#7a7a7a","#cb5c5c","#7a7a7a","#cb5c5c"), 
                    group_palcolor = list(c('#F68282',"#0CB702",'#619CFF')),
                    label_size = 8
) #some aesthetics will be manually modified

png("Output/Re-Batch/600dpi-final/Fig3c-dotplot_overlap_epigene_modifiers_only.png", height = 3.5, width = 5,bg = 'transparent', res = 600, unit = 'in') #4x5.5
ht8$plot
dev.off()

#Figure 3F-G ----
#STEP 1 - obtain DEG list for "stemNK vs. trNK" and "stemNK vs. cytoNK" in B10G5-treated tumors 
##subset B10G5-treated
Idents(NK_fate) <- 'Tx' #levels(NK_fate)
B10G5_subset <- subset(NK_fate, idents = c('B10G5')) #levels(B10G5_subset)

Idents(B10G5_subset) <- 'NK.subset' #levels(B10G5_subset)

stem_tissue <- FindMarkers(B10G5_subset, ident.1 = "stemNK", ident.2 = "trNK", verbose = FALSE, min.pct = 0.25) #default min.pct = 0.01
stem_tissue <- stem_tissue[stem_tissue$p_val_adj <= 0.05,] #need to filter on p.adj <= 0.05 #165
write.csv(stem_tissue,file= "Output/Re-Batch/600dpi-final/03-stem_vs_tissue_DEG.csv",quote=F)

stem_cyto <- FindMarkers(B10G5_subset, ident.1 = "stemNK", ident.2 = "cytoNK", verbose = FALSE,min.pct = 0.25) #default min.pct = 0.01
stem_cyto <- stem_cyto[stem_cyto$p_val_adj <= 0.05,] #need to filter on p.adj <= 0.05 ##84
write.csv(stem_cyto,file= "Output/Re-Batch/600dpi-final/03-stem_vs_cyto_DEG.csv",quote=F)

#STEP 2 - load Kdm6b target list - downloaded from CHIP-Atlas
#(from CHIP-atlas, mm10, Tss+/- 1kb) - compare the two list - 1kb list might be more suitable
Kdm6b_target <- read.delim("MOUSE-Kdm6b.1.tsv",sep="\t") #head(Kdm6b_target) 12297 x 14
Kdm6b_target <- Kdm6b_target[Kdm6b_target$Kdm6b.Average >= 20,]  #set MACS score cut-off - 20 (q value <= 0.01) 11438 x 14
Kdm6b_target <- Kdm6b_target$Target_genes

#STEP 3 - select genes that are both DEG and Kdm6b targets
stem_tissue$gene <- rownames(stem_tissue) 
KDM6b_stem_tissue <- stem_tissue %>%
  filter(gene %in% Kdm6b_target) #dim(KDM6b_stem_tissue)  117 x 6

stem_cyto$gene <- rownames(stem_cyto)
KDM6b_stem_cyto <- stem_cyto %>%
  filter(gene %in% Kdm6b_target) #dim(KDM6b_stem_cyto)  52 x 6

#STEP 4 - volcano plot
library(ggrepel)

##Figure 3F
#stemNK vs trNK
KDM6b_stem_tissue %>% 
  filter(avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) -> top_up

KDM6b_stem_tissue %>% 
  filter(avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 10) -> top_down  

top10 <- rbind(top_up,top_down)  
top10 <- top10$gene

KDM6b_stem_tissue <- KDM6b_stem_tissue %>%
  mutate(text_tf = ifelse(gene %in% top10, gene, NA))

KDM6b_stem_tissue$anno[KDM6b_stem_tissue$avg_log2FC >0 & KDM6b_stem_tissue$gene %in% top10] <- "stemNK"
KDM6b_stem_tissue$anno[KDM6b_stem_tissue$avg_log2FC <0 & KDM6b_stem_tissue$gene %in% top10] <- "trNK"

p1 = ggplot(KDM6b_stem_tissue, aes(x = avg_log2FC, y = -log10(p_val_adj), col = anno)) +
  labs(x = "avg_log2FC", y = "-Log10(P.adj)") +
  theme_classic() +
  geom_vline(xintercept = c(0), col = "gray", linetype = 'dashed') +
  geom_point(alpha = 0.8, shape = 16) + 
  geom_text_repel(aes(label=text_tf), size = 5, max.overlaps = 20) +
  theme(
    axis.title.x = element_text(size = 12),  # Set x-axis label font size and style
    axis.title.y = element_text(size = 12)   # Set y-axis label font size and style
  ) +
  scale_color_manual(values = c('blue','black')) 
#  scale_color_manual(values = c('#0CB702','#619CFF')) 

ggsave(p1,filename = "Output/Re-Batch/600dpi-final/Fig3f-kdm6b_target_stemVtissue.png", device = 'png', height = 4, width = 6, dpi = 600) #Figure 3F

##Figure 3G
#stemNK vs cytoNK
KDM6b_stem_cyto %>% 
  filter(avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) -> top_up

KDM6b_stem_cyto %>% 
  filter(avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 10) -> top_down  

top10 <- rbind(top_up,top_down)  
top10 <- top10$gene
top10 <- c(top10,c('Tcf7'))

KDM6b_stem_cyto <- KDM6b_stem_cyto %>%
  mutate(text_tf = ifelse(gene %in% top10, gene, NA))

KDM6b_stem_cyto$anno[KDM6b_stem_cyto$avg_log2FC >0 & KDM6b_stem_cyto$gene %in% top10] <- "stemNK"
KDM6b_stem_cyto$anno[KDM6b_stem_cyto$avg_log2FC <0 & KDM6b_stem_cyto$gene %in% top10] <- "cytoNK"

p2 = ggplot(KDM6b_stem_cyto, aes(x = avg_log2FC, y = -log10(p_val_adj), col = anno)) +
  labs(x = "avg_log2FC", y = "-Log10(P.adj)") +
  theme_classic() +
  geom_vline(xintercept = c(0), col = "gray", linetype = 'dashed') +
  geom_point(alpha = 0.8, shape = 16) + 
  geom_text_repel(aes(label=text_tf), size = 5, max.overlaps = 20) +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals  
                     limits = c(-6, 6))+
  theme(
    axis.title.x = element_text(size = 12),  # Set x-axis label font size and style
    axis.title.y = element_text(size = 12)   # Set y-axis label font size and style
  ) +
  scale_color_manual(values = c('black','blue')) 
ggsave(p2,filename = "Output/Re-Batch/600dpi-final/Fig3g-kdm6b_target_stemVcyto.png", device = 'png', height = 4, width = 6, dpi = 600) #Figure 3G


#Figure 3J and related- plot human and mouse stemNK overlap ----
length(stemNK.markers) #196 - mouse stemNK marker
TINK <- readRDS(file = "/Users/jem/Documents/Manuscript_Source_Codes/human_TINK/TINK_annotated.rds") #refer to "4-human_TINK_analyses.R"
Idents(TINK) <- 'NK.subset'
hu.markers <- FindAllMarkers(TINK, only.pos = TRUE)
hu.markers <- hu.markers[hu.markers$p_val_adj <= 0.05,]

head(hu.markers) #2222 x 7
hu.stemNK.markers <- hu.markers[hu.markers$cluster == 'stemNK', ] #643 X 7
hu.stemNK.markers <- hu.stemNK.markers$gene #a vector of characters 
length(hu.stemNK.markers) #642

#method 2: use orthology databases 
library(gprofiler2) #https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html

MtoH <- gorth(query = stemNK.markers, source_organism = "mmusculus", 
              target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE,
              numeric_ns = "ENTREZGENE_ACC") #head(MtoH)
M.ortho <- MtoH$ortholog_name #length(M.ortho) #178

overlap3 <- intersect(M.ortho,hu.stemNK.markers)
length(overlap3) #34  

#make venn diagram
y <- list(A = M.ortho,
          B = hu.stemNK.markers)

library(scales) #for color scales

##Figure S7D
plot3 <- ggVennDiagram(y,category.names = c("MOUSE","HUMAN"), label = 'count') + 
  scale_fill_gradient(low = "#effaf6", high = "#a3e1cc") + #set custom color
  theme(legend.position = 'none') + 
  coord_flip() + 
  ggtitle("stemNK markers") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot3, filename = "Output/Re-Batch/600dpi-final/FigS7d-Venn-mouse&human-overlap.png", device = 'png', height = 3, width = 4, dpi = 600) #S7D

#make a dot plot showing only the overlapped chromatin/epigenetic modifiers for human data
markers <- c("SATB1","KDM6B","CHD4")
ph <- DotPlot(TINK, features = markers, dot.scale = 6,cols = c("blue", "red")) + RotatedAxis() + coord_flip()
ggsave(ph, filename = "Output/Re-Batch/600dpi-final/Fig3j-DotPlot_human-mouse-overlapped-epigenetic-modifiers.png", device = 'png', height = 4, width = 5, dpi = 600) #Figure3J
