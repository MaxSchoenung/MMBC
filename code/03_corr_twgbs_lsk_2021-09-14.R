###################################################################
#
# Script: 03_corr_twgbs_lsk_2021-09-14.R
# Correlate the array LSK samples to TWGBS
# Date: 2021-05-12
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
#remotes::install_github("epigen/RnBeads",ref="feature/MouseMethylationBeadChip")
#remotes::install_github("LKremer/ggpointdensity")
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggpointdensity)
library(viridis)
library(ggpubr)
library(methrix)
source("00_utilis.R")

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"

# Load the dataset --------------------------------------------------------
rnb <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name
dim(meth.rnb)

## focus on the cell types
meth.ct <- meth.rnb[,1:27]
pheno.ct <- pheno.rnb[1:27,]

# Load the annotation from stephen ----------------------------------------
annot <- read.delim(paste0(table.dir,"gene-annos_primary_one-row.bed"),stringsAsFactors = F)
annot <- annot %>% 
  rename(
    "Chromosome" = "X.Chromosome")
rownames(annot) <- annot$name
annot.sub <- annot[,c("Chromosome","Start","End","name")]
annot.ranges <- makeGRangesFromDataFrame(annot.sub,keep.extra.columns = T)
seqlevelsStyle(annot.ranges) = "UCSC"
annot.ranges2 <- IRanges::shift(annot.ranges,1) #TWGBS is 1-based and array annotation 0-based


# load twgbs data ---------------------------------------------------------------------
wgbs <- readRDS("/icgc/analysis/OE0565_projects/tce/MouseArray/data/wgbs_mm10_hierachy_methrix/2021-01-29_hierachy_wgbs_all_pops_unfiltered.RDS")
wgbs.lsk <- wgbs[,grep("lsk",colnames(wgbs))[c(1:4,7)]]#use 12347
wgbs.filt <- coverage_filter(m = wgbs.lsk, cov_thr = 10, min_samples = 2) #retained 17,255,601 cpg sites of 20,383,623 (84.65%)
wgbs.mat <- get_matrix(wgbs.filt,"M",add_loci = T,in_granges = T)
wgbs.cov <- get_matrix(wgbs.filt,"C",add_loci = T,in_granges = F)
colMeans(wgbs.cov[,-c(1:3)],na.rm = T)

# Subset by array ---------------------------------------------------------
overlaps <- findOverlaps(wgbs.mat,annot.ranges2,type = "within") #250,007 of 286640 (87.22%)
overlaps.df <- as.data.frame(wgbs.mat[overlaps@from])
overlaps.df$id <- annot.ranges2[overlaps@to]$name
round(nrow(overlaps.df)/nrow(annot)*100,2) #87.22% of all array CpGs
round(nrow(overlaps.df)/nrow(annot[annot$Chromosome%in%paste0("chr",1:19),])*100,2) #93.36% autosomes on mmbc

# Test Plot ---------------------------------------------------------------
wgbs <- data.frame("id"=overlaps.df$id,
                   "wgbs"=rowMeans(overlaps.df[,stringr::str_detect(colnames(overlaps.df),"lsk")],na.rm = T),
                   stringsAsFactors = F)
dim(na.omit(wgbs)) #250,007
dim(wgbs) #250,007
sub.meth <- meth.ct[which(rownames(meth.ct)%in%wgbs$id),]
dim(sub.meth) #234,814
round(nrow(sub.meth)/nrow(annot)*100,2) #81.92% of all array CpGs
round(nrow(sub.meth)/nrow(annot[annot$Chromosome%in%paste0("chr",1:19),])*100,2) #87.69%
array <- data.frame("id"=rownames(sub.meth),
                    "array"=rowMeans(sub.meth[,stringr::str_detect(colnames(sub.meth),"lsk")]),
                    stringsAsFactors = F)
table(wgbs$id%in%array$id)
merged.plot <- merge(wgbs,array,by.x="id",by.y="id",all.x=F,all.y=F)
head(merged.plot)
dim(merged.plot) #234814
table(is.na(merged.plot$wgbs)) #320 NAs

p1 <- print(ggplot(merged.plot,aes(wgbs,array))+
        #geom_pointdensity() +
        stat_pointdensity(geom = GeomPointRast) + 
        scale_color_viridis()+
        theme_classic()+
        ylim(0,1)+
        xlim(0,1)+
        stat_cor()+
        ggtitle("LSK TWGBS"))

#save as pdf
pdf(paste0(plot.dir,Sys.Date(),"_lsk_wgbs_array_corr_raster.pdf"))
p1
dev.off()

# Export the table about coverage information -----------------------------
all.coverage <- get_matrix(wgbs.lsk,"C",add_loci = T,in_granges = T)
overlaps.cov <- findOverlaps(all.coverage,annot.ranges2,type = "within")
overlaps.cov.df <- as.data.frame(all.coverage[overlaps@from,])[,-c(1:5)]


meta.df <- data.frame("sites"=colSums(!is.na(get_matrix(wgbs.lsk,"M",add_loci = T,in_granges = F)[,-c(1:3)])),
                      "median_coverage"=apply(get_matrix(wgbs.lsk,"C",add_loci = T,in_granges = F)[,-c(1:3)],2,function(x)median(x,na.rm = T)),
                      "overlap_array"=colSums(!is.na(overlaps.cov.df)),
                      "median_coverage_array"=apply(overlaps.cov.df,2,function(x)median(x,na.rm = T)),
                      "overlap_cov10"=colSums(overlaps.cov.df>10,na.rm = T),
                      "overlap_cov20"=colSums(overlaps.cov.df>20,na.rm = T)
)
write.table(meta.df,paste0(table.dir,Sys.Date(),"_twgbs_coverage.txt"),sep="\t",quote=F,row.names = T)

## export covered CpGs as a table
write.table(annot[rownames(sub.meth),1:3],paste0(table.dir,Sys.Date(),"_twgbs_overlap_sites.txt"),sep="\t",quote=F,row.names = T)

