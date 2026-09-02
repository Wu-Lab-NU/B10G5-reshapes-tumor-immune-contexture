#KDM6B+H3K27Me3 binding tracks from public datasets (mouse and human)
#For figure S6A-B and figure S7D
#gviz - https://www.bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html
library(Gviz) 
library(GenomicRanges)
library(rtracklayer)  # for importing bigWig and bed

setwd("/Users/jem/Documents/Manuscript_Source_Codes") #update path as needed 

#Load public datasets (mouse)
bw_file <- import.bw("KDM6B+H3K27Me3_binding_inputs/GSM6918459_WT_H3K27me3.bigwig",as = "GRanges") #rtracklayer function- H3K27Me3 in WT mouse NK cells
kdm6b_gr <- import.bed("KDM6B+H3K27Me3_binding_inputs/MOUSE-Oth.ALL.05.Kdm6b.AllCell.bed") #integrated Kdm6b binding data across cell types from CHIP-atlas

#Load public datasets (human)
bw_file <- import.bw("KDM6B+H3K27Me3_binding_inputs/hNK_male_H3K27Me3_ENCFF340WMS.bigWig",as = "GRanges")  #primary
bw_file <- import.bw("KDM6B+H3K27Me3_binding_inputs/hNK_male_H3K27Me3_ENCFF409IKO.bigWig",as = "GRanges")  #back-up - in vitro stimulated NK
kdm6b_gr <- import.bed("KDM6B+H3K27Me3_binding_inputs/HUMAN-Oth.ALL.05.KDM6B.AllCell.bed") 

#Set Genome Region of Interest 
#mouse Tcf7 - chr11:52232682 - 52292036
chrom <- "chr11"
start <- 52232682
end <- 52292036
genome <- "mm10"

#mouse Kdm6b - chr11:69,395,361-69,422,142 
chrom <- "chr11"
start <- 69398361
end <- 69427142 
genome <- "mm10"

#human TCF7 - chr5:134,086,363-134,170,065
chrom <- "chr5"
start <- 134086363
end <- 134170065
genome <- "hg38"

#human TCF7 with enhancers - chr5:134,039,221-134,173,955
chrom <- "chr5"
start <- 134039221
end <- 134173955
genome <- "hg38"

#human KDM6B - chr17:7,816,398-7,854,918
chrom <- "chr17"
start <- 7816398
end <- 7854918
genome <- "hg38"




#H3K27Me3 coverage plot ----
cov_track <- DataTrack(
  range = bw_file,
  genome = genome,
  chromosome = chrom,
  start = start,
  end = end,
  type = "histogram",
  name = " ", #no y axis label for now - leave blank space
  col.histogram = "#2c3e50",
  fill.histogram = "#2c3e50",
  lwd = 0.5, #line width
  ylim = c(0, 5), #pre-check in IGV for range - need to adjust for each species (0-100 for mouse; 0-5 for human)
  aggregation = "mean",     # average values per bin
  window = 100      # bin width in bases (tune this!) - 100 looks good for Tcf7
)

summary(cov_track)


#Kdm6b binding ---- 
kdm6b_sub <- subset(kdm6b_gr, seqnames == chrom & start(kdm6b_gr) < end & end(kdm6b_gr) > start) # Only keep regions in view
kdm6b_merged <- reduce(kdm6b_sub) # Merge overlapping peaks acrosscell types

kdm6b_track <- AnnotationTrack(
  kdm6b_merged,
  genome = genome,
  name = " ", #no y axis label for now
  chromosome = chrom,
  fill = "#2177a4",
  col = "#2177a4",
  cex.title = 0.5
)

#Gene plot with UcscTrack----
gene_track <- UcscTrack(
  genome = genome,
  chromosome = chrom,
  track = "NCBI RefSeq", #"NCBI RefSeq",
  table = "refGene", #glitchy - try ncbiRefSeq or pull out "ucscTables("hg38", "NCBI RefSeq")" (or mm10) to find sub for "refGene" if not working or "ncbiRefSeqCurated"
  from = start,
  to = end,
  trackType = "GeneRegionTrack",
  rstarts = "exonStarts",
  rends = "exonEnds",
  gene = "name2",
  symbol = "name2",
  transcript = "name",
  strand = "strand",
  fill = "#2c3e50",
  col = "#2c3e50",
  name = " ", #no y axis label for now
  cex.title = 1
) 


#(ALTERNATIVE) Gene plot with GeneRegionTrack----
##note this is a temp fix that replcaes the original UcscTrack() lines, which are not working as of Sep 2026....
q <- ucscTableQuery(
  "hg38",
  table = "ncbiRefSeq",
  range = GRanges(
    chrom,
    IRanges(start, end)
  )
)

refseq_test <- getTable(q)

gene_track <- GeneRegionTrack(
  chromosome = chrom,
  genome = "hg38",
  start = start,
  end = end,
  rstarts = refseq_test$exonStarts,
  rends = refseq_test$exonEnds,
  gene = refseq_test$name2,
  symbol = refseq_test$name2,
  transcript = refseq_test$name,
  strand = refseq_test$strand,
  fill = "#2c3e50",
  col = "#2c3e50",
  name = " ",
  cex.title = 1
)

#(ALTERNATIVE) gene plot forward strand only ---
refseq_fwd <- refseq_test[refseq_test$strand == "+", ]

gene_track_fwd <- GeneRegionTrack(
  chromosome = chrom,
  genome = "hg38",
  start = start,
  end = end,
  rstarts = refseq_fwd$exonStarts,
  rends = refseq_fwd$exonEnds,
  gene = refseq_fwd$name2,
  symbol = refseq_fwd$name2,
  transcript = refseq_fwd$name,
  strand = refseq_fwd$strand,
  
  collapseTranscripts = "meta",
  stacking = "squish",
  showId = TRUE,
  transcriptAnnotation = "symbol",
  
  fill = "#2c3e50",
  col = "#2c3e50",
  name = " "
)

#set visuals for gene_track
displayPars(gene_track) <- list(
  fill = "#2c3e50",          # deep blue fill for exons
  col = NA,                  # no exon borders
  col.line = "gray40",       # intron arrow line color
  collapseTranscripts = "meta",
  transcriptAnnotation = "symbol",
  showId = T,
  just.group = "below",      # place symbols above transcripts
  cex.group = 1,                 # label font size
  fontface = 1,              # plain font
  arrowHeadWidth = 8         # make arrows cleaner
)

plotTracks(gene_track) #check as sometimes glitchy

#highlight KDM6B binding regions ----
highlight_track <- HighlightTrack(
  trackList = list(cov_track),  # tracks to overlay highlights on
  start = start(kdm6b_merged),
  end = end(kdm6b_merged),
  chromosome = chrom,
  genome = genome,
  col = NA,
  fill = "#2177a460"  
)


#genomic range axis ----
axis_track <- GenomeAxisTrack(
  genome = genome,
  chromosome = chrom,
  name = "",                 # no label
  col = "black",             # axis line color
  col.line = "black",        # tick marks
  fontcolor = "black",       # label color
  cex = 0.8,                 # font size
  littleTicks = F,        # small ticks between major ones
  lwd = 0.5,
  labelPos = "alternating",         # place labels below the axis line
  exponent = 6 #3 = kb,  6 = mb
)


#test graph before finalization ----
plotTracks(
  list(highlight_track, kdm6b_track, gene_track_fwd,axis_track), #gene_track_fwd
  from = start,
  to = end,
  cex.title = 0.8,
  cex.axis = 0.8,
  col.axis = "black",
  type = 'histogram'
) 

#assemble, plot and save ---- 
#change file name/gene name as needed
png("Output/Fig3/FigS6b-supp-human-tcf7-track+enhancers.png",width = 4.5, height = 3, res = 600, unit = 'in')  # 5x3 adjust dimensions - may see error message if too small to fit the graph

plotTracks(
  list(highlight_track, kdm6b_track, gene_track_fwd,axis_track),
  from = start,
  to = end,
  transcriptAnnotation = "symbol",
  collapseTranscripts = "meta",  # ← THIS collapses isoforms into one
  background.title = "white",
  col.title = "black",
  fontcolor.title = "black",
  cex.title = 0.8,
  cex.axis = 0.8,
  col.axis = "black",
  type = 'histogram',
  sizes = c(3,0.5,0.5,0.8)
)

dev.off()


