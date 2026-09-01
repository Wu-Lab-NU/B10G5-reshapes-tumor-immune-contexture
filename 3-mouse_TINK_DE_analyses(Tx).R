#DE analysis comparing B10G5-treated and untreated tumor-infiltrating NK cells 
#EnrichR data used for Figure 2E
#B10G5 up DEG used for Figure 3A

#PREP - skip this section if continued from "1-mouse_TINK_basic_analysis.R
library(Seurat)
library(dplyr)

setwd("/Users/jem/Documents/Manuscript_Source_Codes") #update path as needed 
NK_fate <- readRDS("NK_data_batch.rds") #read NK object

# find markers between Tx and B10G5-treated
Idents(NK_fate) <- 'Tx'
markers1 <- FindMarkers(NK_fate, ident.1 = "B10G5", ident.2 = "No Tx", verbose = FALSE)
markers1 <- markers1[markers1$p_val_adj <= 0.05,] #need to filter on p.adj <= 0.05

#pathway enrichment analysis 

#enrichR
#perform EnrichR pathway enrichment analysis 
#https://cran.r-project.org/web/packages/enrichR/readme/README.html
library(enrichR)

#connect to the online database
#By default human genes are selected otherwise select your organism of choice.
listEnrichrSites()

#find the list of all available databases from Enrichr.
dbs <- listEnrichrDbs()
head(dbs)

#View and select your favourite databases
dbs <- c("KEGG_2019_Mouse","Reactome_2022","MSigDB_Hallmark_2020","ENCODE_TF_ChIP-seq_2015","GO_Biological_Process_2025","GO_Biological_Process_2015")

#import data into data frames
DEG = markers1

#order based on log2Fc, from small to big 
DEG = DEG[order(DEG$avg_log2FC),]

#get the up and down regulated DEG list (up = up in B10G5 treated)
up = rownames(DEG[DEG$avg_log2FC > 0,])
down = rownames(DEG[DEG$avg_log2FC < 0, ])

#perform enrichR --> output: a large list containing dfs for each database
up.enrich = enrichr(up, dbs)
# Wait for 5 s to avoid overloading the server
Sys.sleep(5)
#!!!! the server will return a wrong list if you do not wait here
down.enrich = enrichr(down, dbs)
# Wait for 5 s to avoid overloading the server
Sys.sleep(5)

head(up.enrich[[1]])
head(down.enrich[[1]])

#export the tables in formats ready for "PathwayBarGraph.R"
#save up- and down- pathways in one csv file + column "Categories" specifying UP vs DOWN 
for (j in 1:length(dbs)) {
  #add categories UP or DOWN
  up.term = up.enrich[[j]]
  down.term = down.enrich[[j]]
  up.term$Categories = "UP"
  down.term$Categories = "DOWN"
  
  #remove empty columns 
  up.term = subset(up.term, select = -c(Old.P.value,Old.Adjusted.P.value))
  down.term = subset(down.term, select = -c(Old.P.value,Old.Adjusted.P.value))
  
  #filter terms with p.adj <= 0.05
  up.term = subset(up.term, Adjusted.P.value <= 0.05)
  down.term = subset(down.term, Adjusted.P.value <= 0.05)
  
  #combine the two dfs and export
  enrich <- rbind(up.term,down.term)
  filename = paste("Output/Re-Batch/EnrichR/B10G5_v_noTx_",dbs[j],".xlsx", sep="")
  write.xlsx(enrich,file= filename,rowNames = T) 
  
}
##end of the entire for loop####
