###################################################################
#
# Script: 05_corr_rrbs.R
# Correlation to RRBS data of C. Bock (PMID: 22841485)
# Date: 2021-05-13
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
library(gridExtra)


# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Set the colours ---------------------------------------------------------
colours <- c("LSK"="gray50",
             "CMP"="darkorange",
             "GMP"="darkgoldenrod2",
             "MEP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "CD8_t-cell"="lightsteelblue1",
             "CD4_t-cell"="lightskyblue")

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

## make GRanges
annot.sub <- annot[,c("Chromosome","Start","End","name")]
annot.ranges <- makeGRangesFromDataFrame(annot.sub,keep.extra.columns = T)
seqlevelsStyle(annot.ranges) = "UCSC" 

# load rrbs data ---------------------------------------------------------------------
## function for shaping the files
shape_fun <- function(x){
  v4_rm <- stringr::str_remove_all(x$V4,"'")
  v4_m <- stringr::str_split_fixed(v4_rm,"/",2)[,1]
  v4_u <- stringr::str_split_fixed(v4_rm,"/",2)[,2]
  return.df <- data.frame("ind"=paste0(x$V1,"_",x$V2,"_",x$V3),
                          "chr"=x$V1,
                          "start"=x$V2,
                          "end"=x$V3,
                          "strand"=x$V6,
                          "U"=v4_u,
                          "M"=v4_m,
                          "cov"=(as.numeric(v4_u)+as.numeric(v4_m)),
                          "beta"=x$V5)
  
}


files.load <- paste0(data.dir,"bock_rrbs_2012/",list.files(paste0(data.dir,"bock_rrbs_2012/")))
index <- stringr::str_remove_all(list.files(paste0(data.dir,"bock_rrbs_2012/")),"RRBS_cpgMethylation_")
index <- stringr::str_remove_all(index,".bed")
file.list <- list()
shaped.list <- list()
for(i in 1:length(files.load)){
  message(paste0("Processing: ",index[i]))
  message(paste0("File: ",files.load[i]))
  file.list[[index[i]]] <- read.delim(files.load[i],header = F)
  shaped.list[[index[i]]] <- shape_fun(file.list[[index[i]]])
}

## document the meta data for the rrbs data
meta.df <- data.frame("sites"=do.call(c,lapply(file.list,nrow)))

# autosomes filter
list.filt <- lapply(shaped.list,function(x){x[x$chr%in%paste0("chr",1:19),]})
meta.df$n_autosomes <- do.call(c,lapply(list.filt,nrow))

# Analyze just those sites with a 20x coverage ----------------------------
list.filt <- lapply(list.filt,function(x){x[x$cov>20,]})
meta.df$cov_passed <- do.call(c,lapply(list.filt,nrow))

## merge the beta values of the list
beta.df <- list.filt[[1]][,c("ind","beta")]
colnames(beta.df) <- c("ind",index[1])
for(i in 2:length(list.filt)){
  message("Processing: ",names(list.filt)[i])
  beta.df <- merge(beta.df,list.filt[[i]][,c("ind","beta")],by.x="ind",by.y="ind",all.x=T,all.y=T,sort=F)
  colnames(beta.df)[i+1] <- names(list.filt)[i]
}
(nrow(beta.df)-colSums(is.na(beta.df)))[-1]==meta.df$cov_passed #check!


# Make GRanges ---------------------------------------------------------------
beta.sub <- beta.df[,-1]
beta.sub$chr <- stringr::str_split_fixed(beta.df[,1],pattern = "_",3)[,1]
beta.sub$start <- stringr::str_split_fixed(beta.df[,1],pattern = "_",3)[,2]
beta.sub$end <- stringr::str_split_fixed(beta.df[,1],pattern = "_",3)[,3]
beta.ranges <- makeGRangesFromDataFrame(beta.sub,keep.extra.columns = T)


# Liftover ----------------------------------------------------------------
ch <- import.chain(paste0(table.dir,"mm9ToMm10.over.chain"))
seqlevelsStyle(beta.ranges) = "UCSC" 
beta.mm10 <-  liftOver(beta.ranges, ch)
beta.mm10 <- unlist(beta.mm10)
genome(beta.mm10) = "mm10"
width(beta.mm10) <- 2 #focus on C as in annotation


# subset with annot -------------------------------------------------------
overlaps <- findOverlaps(beta.mm10,annot.ranges,type = "within")
overlaps.df <- as.data.frame(beta.mm10[overlaps@from])
overlaps.df$id <- annot.ranges[overlaps@to]$name

meta.df$on_array <- (nrow(overlaps.df)-colSums(is.na(overlaps.df)))[-c(1:5,22)]


rrbs.names <- c("CMP","GMP","MEP","CD4","CD8","B_cell","Gran","Mono")
array.names <- c("cmp","gmp","mep","cd4","cd8","bcell","neutro","mono")

# Plot ---------------------------------------------------------------
plot.list <- list()

for(i in 1:length(rrbs.names)){
rrbs <- data.frame("id"=overlaps.df$id,
                   "rrbs"=rowMeans(overlaps.df[,stringr::str_detect(colnames(overlaps.df),rrbs.names[i])])/1000,
                   stringsAsFactors = F)
sub.meth <- meth.ct[which(rownames(meth.ct)%in%rrbs$id),]
array <- data.frame("id"=rownames(sub.meth),
                   "array"=rowMeans(sub.meth[,stringr::str_detect(colnames(sub.meth),array.names[i])]),
                   stringsAsFactors = F)

merged.plot <- merge(rrbs,array,by.x="id",by.y="id",all.x=F,all.y=F)
merged.plot$ct <- rep(array.names[i],nrow(merged.plot))
plot.list[[i]] <- merged.plot
}

plot.df <- do.call(rbind,plot.list)

pdf(paste0(plot.dir,Sys.Date(),"_rrbs_array_corr_raster_grid.pdf"),width = 12,height = 10)
ggplot(plot.df,aes(rrbs,array))+
  stat_pointdensity(geom = GeomPointRast) + 
  facet_wrap(~ct)+
  scale_color_viridis()+
  theme_classic()+
  ylim(0,1)+
  xlim(0,1)+
  stat_cor()+
  geom_text(data=plyr::count(na.omit(plot.df), vars = c("ct")),
            aes(x=0.1, y=0.8, label=paste0("n=",freq)),
            colour="black", inherit.aes=FALSE, parse=FALSE)
dev.off()