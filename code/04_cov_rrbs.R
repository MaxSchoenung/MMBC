###################################################################
#
# Script: 04_cov_rrbs.R
# Coverage of Array CpGs in RRBS data
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
#remotes::install_github("epigen/RnBeads",ref="feature/MouseMethylationBeadChip")
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)
library(dplyr)
library(purrr)
library(GenomicRanges)
library(rtracklayer)
library(gridExtra)
library(grid)
library(gridtext)

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

# Load the annotation from stephen ----------------------------------------
annot <- read.delim(paste0(table.dir,"gene-annos_primary_one-row.bed"),stringsAsFactors = F)
annot <- annot %>% 
  rename(
    "Chromosome" = "X.Chromosome")
rownames(annot) <- annot$name

## convert to GRanges
annot.sub <- annot[,c("Chromosome","Start","End","name")]
annot.ranges <- makeGRangesFromDataFrame(annot.sub,keep.extra.columns = T)
seqlevelsStyle(annot.ranges) = "UCSC" 

## import liftOver chain
ch <- import.chain(paste0(table.dir,"mm9ToMm10.over.chain"))

# load rrbs data ---------------------------------------------------------------------
## function for shaping the files
shape_fun <- function(x){
  v4_rm <- stringr::str_remove_all(x$V4,"'")
  v4_m <- stringr::str_split_fixed(v4_rm,"/",2)[,1]
  v4_u <- stringr::str_split_fixed(v4_rm,"/",2)[,2]
  return.df <- data.frame("chr"=x$V1,
                          "start"=x$V2,
                          "end"=x$V3,
                          "strand"=x$V6,
                          "U"=v4_u,
                          "M"=v4_m,
                          "cov"=(as.numeric(v4_u)+as.numeric(v4_m)),
                          "beta"=x$V5)
  return.df <- return.df[return.df$chr%in%paste0("chr",1:19),] #autosomes filter!!!
  message(paste0(nrow(return.df)," sites detected."))
  message(paste0(sum(return.df$cov>20)," sites with coverage > 20."))
  return(return.df)
}

## function for GRanges and liftOver
lift_fun <- function(x){
  ch <- import.chain(paste0(table.dir,"mm9ToMm10.over.chain"))
  beta.ranges <- makeGRangesFromDataFrame(x,keep.extra.columns = T)
  width(beta.ranges) <- 2 ##bock data is including CG but we wanna focus on C
  seqlevelsStyle(beta.ranges) = "UCSC" 
  beta.mm10 <-  liftOver(beta.ranges, ch)
  beta.mm10 <- unlist(beta.mm10)
  genome(beta.mm10) = "mm10"
  return(beta.mm10)
}

## overlaps fun
over_fun <- function(x){
  overlaps <- findOverlaps(x,annot.ranges,type = "within")
  over.df <- x[overlaps@from]
  over.df$probe <- annot.sub[overlaps@to,"name"]
  return(over.df)
}



# Process the RRBS data ---------------------------------------------------
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

## liftover
mm10.list <- lapply(shaped.list, lift_fun)

## overlaps
over.list <- lapply(mm10.list,over_fun)

## extract coverage
cov.list <- lapply(over.list,function(x)x$cov)
lapply(cov.list,median)

plot.list <- list()
## plotting
for(i in 1:length(cov.list)){
  plot.df <- data.frame("value"=cov.list[[i]])
  plot.list[[i]] <- ggplot(plot.df,aes(x=value))+
    geom_histogram(bins=100)+
    xlim(0,200)+
    theme_classic()+
    xlab("Coverage")+
    ggtitle(paste0(names(cov.list)[i]," (n=",nrow(na.omit(plot.df))," sites)"))+
    theme(plot.title = element_text(size = 10, face = "bold"))
}

# add common axis labels
plot.list <- plot.list%>% map(~.x + labs(x=NULL, y=NULL))

# plotmath expressions
yleft <- textGrob("Count", 
                  rot = 90, gp = gpar(fontsize = 15,fontface=2))

bottom <- textGrob("Coverage", gp = gpar(fontsize = 15,fontface=2))


pdf(paste0(plot.dir,Sys.Date(),"_hist_rrbs_cov_list.pdf"),width = 9,height = 9)
grid.arrange(grobs = plot.list, ncol = 4,
             left=yleft,bottom=bottom)
dev.off()

## coverage meta data
meta.df <- data.frame("sites"=do.call(c,lapply(shaped.list,nrow)),
                      "median_coverage"=do.call(c,lapply(shaped.list,function(x)median(x$cov))),
                      "overlap_array"=do.call(c,lapply(over.list,length)),
                      "median_coverage_array"=do.call(c,lapply(over.list,function(x)median(x$cov))),
                      "overlap_cov10"=do.call(c,lapply(over.list,function(x)sum(x$cov>10))),
                      "overlap_cov20"=do.call(c,lapply(over.list,function(x)sum(x$cov>20)))
                      )
write.table(meta.df,paste0(table.dir,Sys.Date(),"_rrbs_coverage.txt"),sep="\t",quote=F,row.names = T)

## median cov all overlapping samples
do.call(c,lapply(over.list,function(x)x$cov))%>%median
        