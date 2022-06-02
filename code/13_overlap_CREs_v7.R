###################################################################
#
# Script: 13_overlap_CREs_v7.R
# Overlap between CRE catalogues and diffDMPs
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
library(dplyr)
library(pheatmap)
library(biomaRt)
library(RnBeads)
library(GenomicRanges)
library(ggplot2)
library(VennDiagram)
library(UpSetR)
library(ComplexHeatmap)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Define the colours ------------------------------------------------------
colours <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "CD8_t-cell"="lightsteelblue1",
             "CD4_t-cell"="lightskyblue")

# Load the methylation data -----------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)
rm(rnb.full)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name

# Load the HSC dmps -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
dmps.unique <- dmps[!duplicated(dmps$site),] #37,512 unique sites

# Load the annotation -----------------------------------------------------
annot.full <- read.delim(paste0(table.dir,"2022-05-25_MMBC_annotation_modified.bed"),stringsAsFactors = F)
rownames(annot.full) <- annot.full$IlmnID
annot.full <- annot.full[annot.full$Chromosome%in%paste0("chr",1:19),]
annot <- annot.full[,1:3]
annot.ranges <- GenomicRanges::makeGRangesFromDataFrame(annot)
seqlevelsStyle(annot.ranges) = "UCSC"

## annotate the dmps
dmps.annot <- cbind(annot[dmps.unique$site,],dmps.unique)
dmp.ranges <- makeGRangesFromDataFrame(dmps.annot,keep.extra.columns = T)

# Shape the enhancer ------------------------------------------------------
enhancer <- read.delim(paste0(table.dir,"amit2014-enhancers-v1_cluster-all_mm10.bed"),stringsAsFactors = F)
colnames(enhancer)[1] <- "Chromosome"
enhancer <- enhancer[which(enhancer$Chromosome%in%c(1:19)),]

enhancer$ClusterID[enhancer$ClusterID == 1] <- "common"
enhancer$ClusterID[enhancer$ClusterID == 2] <- "lymphoid_progenitors"
enhancer$ClusterID[enhancer$ClusterID == 3] <- "T_NK_cells"
enhancer$ClusterID[enhancer$ClusterID == 4] <- "B_cells"
enhancer$ClusterID[enhancer$ClusterID == 5] <- "erythroid"
enhancer$ClusterID[enhancer$ClusterID == 6] <- "myeloid_progenitors"
enhancer$ClusterID[enhancer$ClusterID == 7] <- "myeloid"
enhancer$ClusterID[enhancer$ClusterID == 8] <- "progenitors"
enhancer$ClusterID[enhancer$ClusterID == 9] <- "erythroid_progenitors"


enhancer.ranges <- GenomicRanges::makeGRangesFromDataFrame(enhancer, seqnames.field = "Chromosome",start.field = "Start",end.field = "End",keep.extra.columns = T)
seqlevelsStyle(enhancer.ranges) <- "UCSC"
unique(enhancer.ranges)

# Immune CRE --------------------------------------------------------------
immune.cre <- read.csv(paste0(table.dir,"ImmGenATAC18_AllOCRsInfo.csv"))
immune.cre.auto <- immune.cre[immune.cre$chrom%in%paste0("chr",1:19),]
immune.cre.ranges <- GenomicRanges::makeGRangesFromDataFrame(immune.cre.auto,seqnames.field = "chrom",start.field="Summit",end.field="Summit")
immune.cre.ranges.ext <- immune.cre.ranges+125 #as the summits are 250bp centered
seqlevelsStyle(immune.cre.ranges.ext) <- "UCSC"

# Load the cCRE --------------------------------------------------------
ccre <- read.delim(paste0(table.dir,"mm10_cCRE"))
ccre2 <- ccre[,c("X.chrom","chromStart","chromEnd","encodeLabel")]
ccre2 <- ccre2[ccre2$X.chrom%in%paste0("chr",1:19),]
ccre.ranges <- makeGRangesFromDataFrame(ccre2, seqnames.field = "X.chrom",start.field = "chromStart",end.field = "chromEnd",keep.extra.columns = T,ignore.strand = T)
seqlevelsStyle(ccre.ranges) <- "UCSC"

# VISION CRE/gene pairs -----------------------------------------------------
# info: https://ccre-filter.bx.psu.edu/Notes.txt
#vision.cre.gene <- read.delim("http://usevision.org/data/mm10/cCRE-gene-pairs/mouse_allChr_wTAD_weRP_wAllComp_071819.txt") #CRE gene pairs
vision.cre <- read.delim("https://usevision.org/data/mm10/VISIONmusHem_ccREs.txt",header = F) #CRE
vision.cre2 <- vision.cre[vision.cre$V1%in%paste0("chr",1:19),]
vision.cre.ranges <- makeGRangesFromDataFrame(vision.cre2, seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T,ignore.strand = T)
seqlevelsStyle(vision.cre.ranges) <- "UCSC"

# Just use probes on array ------------------------------------------------
peak.list <- list("cCRE"=ccre.ranges,"Enhancer"=enhancer.ranges,"ImmuneCRE"=immune.cre.ranges.ext,"VisionCRE"=vision.cre.ranges)
array.peaks <- lapply(peak.list,function(x)findOverlaps(annot.ranges,x))
plot.over.df <- cbind("Overlap"=sapply(array.peaks,function(x)length(unique(x@to))),
      "Total"=sapply(array.peaks,function(x)x@nRnode))%>%as.data.frame()
plot.over.df$Unique <- plot.over.df$Total-plot.over.df$Overlap
plot.over.df.melted <- reshape2::melt(as.matrix(plot.over.df[,c(1,3)]))
pdf(paste0(plot.dir,Sys.Date(),"_mmbc_probes_overlap_array.pdf"))
ggplot(plot.over.df.melted,aes(Var1,value,fill=Var2))+
  geom_bar(stat="identity",colour="black")+
  theme_classic()+
  ylab("CREs [n]")+
  xlab("Catalogue")+
  scale_fill_grey()+
  scale_y_continuous(labels = scales::comma)
dev.off()

# Overlap between catalogues and DMPs -------------------------------------
peak.list.all <- list("DMPs"=dmp.ranges,"cCRE"=ccre.ranges,"Enhancer"=enhancer.ranges,"ImmuneCRE"=immune.cre.ranges.ext,"VisionCRE"=vision.cre.ranges)
peak.names <- lapply(peak.list.all,function(x)unique(names(annot.ranges)[findOverlaps(annot.ranges,x)@from]))
peak.names$DMPs <- peak.names$DMPs[peak.names$DMPs%in%names(dmp.ranges)] #some of the DMPs have multiple overlaps
lapply(peak.names,length)

## prefered: upset plot
pdf(paste0(plot.dir,Sys.Date(),"_CREs_upset.pdf"))
upset(fromList(peak.names), order.by = "freq")
dev.off()

## extract peaks of chromCRE, previous mdCRE and novel mdCRE
previousCRE <- do.call(c,peak.names[-1])
chromCRE <- unique(previousCRE[!previousCRE%in%peak.names$DMPs])
novel_mdCRE <- unique(peak.names$DMPs[!peak.names$DMPs%in%previousCRE])
overlap_mdCRE <- unique(peak.names$DMPs[peak.names$DMPs%in%previousCRE])

length(chromCRE) #94,943
length(novel_mdCRE) #12,856
length(overlap_mdCRE) #24,656


# Check -------------------------------------------------------------------
mdCRE2 <- readRDS(paste0(data.dir,"2022-03-08_overlapping_CRE.RDS"))
chromCRE2 <- readRDS(paste0(data.dir,"2022-03-08_previous_CRE.RDS"))
novel_mdCRE2 <- readRDS(paste0(data.dir,"2022-03-08_novel_mdCRE.RDS"))
all(mdCRE2%in%overlap_mdCRE)
all(chromCRE2%in%chromCRE)
all(novel_mdCRE2%in%novel_mdCRE)

pdf(paste0(plot.dir,Sys.Date(),"_CRE_Venn.pdf"))
plot(eulerr::euler(c(A=length(chromCRE), B=length(novel_mdCRE), "A&B"=length(overlap_mdCRE))))
upset(fromList(list("previousCRE"=previousCRE,"DMPs"=peak.names$DMPs)), order.by = "freq")
dev.off()

# Save the sets -----------------------------------------------------------
saveRDS(novel_mdCRE,paste0(data.dir,Sys.Date(),"_novel_mdCRE.RDS"))
saveRDS(overlap_mdCRE,paste0(data.dir,Sys.Date(),"_overlapping_CRE.RDS"))
saveRDS(chromCRE,paste0(data.dir,Sys.Date(),"_previous_CRE.RDS"))

# Add these sets to peak names --------------------------------------------
peak.names.ext <- c(peak.names,list("previous_mdCRE"=overlap_mdCRE,"chromCRE"=chromCRE,"cmdCRE"=novel_mdCRE,"all_array"=annot.full$name))

# Analyze the feature class distribution of those CREs --------------------
feat.class <- lapply(peak.names.ext,function(x)round(table(annot.full[x,"Genomic_region_class"])/length(x)*100,2))
feat.class <- do.call(rbind,feat.class)[,c(3,2,7,4,6,1,5)]
pdf(paste0(plot.dir,Sys.Date(),"_CREs_feat_class.pdf"))
pheatmap(feat.class,
         breaks=seq(0, 50, by = 1),
         cluster_rows = T,
         cluster_cols = F,
         treeheight_row = 0,
         display_numbers = T,
         cellheight = 16,
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 50, by = 1))))
pheatmap(feat.class[1:5,],
         breaks=seq(0, 50, by = 1),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         cellheight = 16,
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 50, by = 1))))
pheatmap(feat.class[6:8,],
         breaks=seq(0, 50, by = 1),
         cluster_rows = T,
         cluster_cols = F,
         display_numbers = T,
         treeheight_row = 0,
         cellheight = 16,
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 50, by = 1))))
dev.off()


# Analyze the distribution over regulatory built --------------------------
# Read the regulatory build
reg.build <- ape::read.gff(paste0(table.dir,"mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff"))
reg.build$feat_class <- stringr::str_split_fixed(reg.build$attributes,"feature_type=",2)[,2]
reg.build <- reg.build[reg.build$seqid%in%c(1:19),]
reg.ranges <- GenomicRanges::makeGRangesFromDataFrame(reg.build, seqnames.field = "seqid",start.field = "start",end.field = "end",keep.extra.columns = T)
seqlevelsStyle(reg.ranges) <- "UCSC"
reg.ranges <- sort(reg.ranges)

# Add to annotation
over <- findOverlaps(reg.ranges,annot.ranges)
over.df <- as.data.frame(over)
over.df$class <- reg.ranges$feat_class[over.df$queryHits]
over.df$site <- names(annot.ranges)[over.df$subjectHits]

# Plot the classes
feat.class <- lapply(peak.names.ext,function(x)round(table(over.df[which(over.df$site%in%x),"class"])/length(x)*100,2))
feat.class <- do.call(rbind,feat.class)
pdf(paste0(plot.dir,Sys.Date(),"_CREs_reg_built.pdf"))
pheatmap(feat.class,
         breaks=seq(0, 50, by = 1),
         cluster_rows = T,
         cluster_cols = F,
         treeheight_row = 0,
         display_numbers = T,
         cellheight = 16,
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 50, by = 1))))
pheatmap(feat.class[1:5,],
         breaks=seq(0, 50, by = 1),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         cellheight = 16,
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 50, by = 1))))
pheatmap(feat.class[6:8,],
         breaks=seq(0, 50, by = 1),
         cluster_rows = T,
         cluster_cols = F,
         display_numbers = T,
         treeheight_row = 0,
         cellheight = 16,
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 50, by = 1))))
dev.off()


# Overlap with MMBC ------------------------------------
## generate the annotation
colours <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "CD8_t-cell"="lightsteelblue1",
             "CD4_t-cell"="lightskyblue")

anno.df <- pheno.rnb[,"CellType",drop=F]
rownames(anno.df) <- pheno.rnb[,"Sample_Name"]

order.vec <- c(4:6,1:3,7:12,22:27,19:21,13:18)
pdf(paste0(plot.dir,Sys.Date(),"_cCRE_corr.pdf"))
for(i in 1:length(peak.names.ext)){
corr.df <- meth.rnb[rownames(meth.rnb)%in%peak.names.ext[[i]],]%>%
    cor
print(
  pheatmap(
    corr.df[order.vec,order.vec],
         color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0.8, 1, by = 0.001))),
         breaks=seq(0.8, 1, by = 0.001),
         cluster_rows = F,
         cluster_cols = F,
         border_color = NA,
         show_colnames = F,
         show_rownames=F,
         annotation_col = anno.df,
         annotation_row = anno.df,
         annotation_colors = list("CellType"=colours),
         main=names(peak.names.ext)[i]))
print(
  pheatmap(
    corr.df[order.vec,order.vec],
    color = colorRampPalette(c("white", "lavenderblush","darkred"))(length(seq(0, 1, by = 0.01))),
    breaks=seq(0, 1, by = 0.01),
    cluster_rows = F,
    cluster_cols = F,
    border_color = NA,
    main=names(peak.names.ext)[i]))
}
dev.off()


