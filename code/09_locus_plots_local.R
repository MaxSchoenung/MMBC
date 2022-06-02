###################################################################
#
# Script: 09_locus_plots_local.R
# Locus plots for cell type marker
# Date: 2021-07-12
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
library(dplyr)
library(Gviz)
library(biomaRt)
library(GenomicRanges)

# Set Directories ---------------------------------------------------------
setwd("~/Desktop/")
table.dir <- "/Users/maximilianschoenung/Documents/AG_Lipka/Projects/Mouse_CellType_Array/Tables/"
plot.dir <- "/Users/maximilianschoenung/Documents/AG_Lipka/Projects/Mouse_CellType_Array/Plots/"
data.dir <- "/Users/maximilianschoenung/Documents/AG_Lipka/Projects/Mouse_CellType_Array/Data/"


# Define the biomart ------------------------------------------------------
mm10 <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="http://nov2020.archive.ensembl.org") #biomart for annotation

# Load the HSC dmps -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
dmps[dmps$site=="cg45931459_TC21",]

# Load the methylation data -----------------------------------------------
mean.meth <- readRDS(paste0(data.dir,"2021-07-13_mean_meth_allPops.RDS"))
mean.meth <- mean.meth[,c(1:4,8,9,5:7)]

# Load the annotation from stephen ----------------------------------------
annot <- read.delim(paste0(table.dir,"Annotations/2022-05-25_MMBC_annotation_modified.bed"),stringsAsFactors = F)
rownames(annot) <- annot$IlmnID
annot.sub <- annot[annot$rnbeads==T,c("Chromosome","Start","End","IlmnID")]
annot.ranges <- GenomicRanges::makeGRangesFromDataFrame(annot.sub,keep.extra.columns = T)
seqlevelsStyle(annot.ranges) = "UCSC"
unique.ranges <- annot.ranges[unique(dmps$site),]

# Define the colours ------------------------------------------------------
colours <- c("LSK"="gray50",
             "MEP"="orangered4",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "CD4"="lightsteelblue1",
             "CD8"="lightskyblue",
             "B-cell"="navy")

################################
#### Select Gene and Range #####
# ################################
selected.symbol="Cebpe"
down <- 2000
up <- 1000
# ###############################


## Download the TSS for gene to define the region
tx.id <- getBM(attributes=c("ensembl_gene_id","external_gene_name","ensembl_transcript_id_version","chromosome_name","transcription_start_site"),filters="external_gene_name",values=selected.symbol,mart=mm10)
tx.id <- tx.id[tx.id$chromosome_name%in%1:19,]
tx.ranges <- makeGRangesFromDataFrame(tx.id, seqnames.field = "chromosome_name",start.field = "transcription_start_site",end.field = "transcription_start_site",keep.extra.columns = T,ignore.strand = T)
seqlevelsStyle(tx.ranges) <- "UCSC"

selected.region <- tx.ranges[tx.ranges$external_gene_name==selected.symbol,][1]
chr <- as.character(unique(seqnames(selected.region)))
gen <- "mm10"
from <- selected.region@ranges@start-down
to <- selected.region@ranges@start+up

# Generate the tracks -----------------------------------------------------
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
gtrack <- GenomeAxisTrack(labelPos = "above")
dmptrack <- AnnotationTrack(unique.ranges,chromosome = chr, name = "DMPs")

# add the methylation data
overl <- findOverlaps(annot.ranges,GRanges(paste0(chr,":",from,"-",to)))
meth.tracks <- list()
for(i in 1:ncol(mean.meth)){
  track <- annot.ranges[overl@from,]
  mcols(track)$meth <- mean.meth[names(track),i]
  meth.tracks[[i]] <- DataTrack(track, name = colnames(mean.meth)[i],type="h",col=colours[i],
                                ylim=c(0,1.01))
}

biomTrack <- BiomartGeneRegionTrack(genome = "mm10", name = "ENSEMBL", 
                                    symbol = selected.symbol, biomart = mm10)
grtrack <- GeneRegionTrack(biomTrack@range, genome = gen,
                           chromosome = chr, name = "Gene Model")
## change visualization pars
plots <- c(list(itrack,gtrack,grtrack),meth.tracks,list(dmptrack))

pdf(paste0(plot.dir,Sys.Date(),"_gviz_cebpe.pdf"),height = 5)
plotTracks(plots,
           from = from, to = to,transcriptAnnotation = "symbol",#sizes=c(1,1,1,rep(1,9),1),
           background.panel = "white", background.title = "white",
           col.axis="black",fontcolor.title="black")
dev.off()