###################################################################
#
# Script: 01_preprocessing.R
# Import of Mouse Methylation Array Data into RnBeads
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

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
report.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/reports/2021/"
sample_sheet <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/sample_sheet_schoenun_17.02.2021.txt"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"

# RnBeads Options ---------------------------------------------------------
rnb.options(
  identifiers.column = "Sample_Name",
  normalization.method = "scaling",
  normalization.background.method = "subtraction",
  filtering.greedycut.pvalue.threshold = 0.01,
  filtering.sex.chromosomes.removal = TRUE, 
  filtering.missing.value.quantile = 0, #remove all NAs
  disk.dump.big.matrices = FALSE,
  import.table.separator="\t"
)

# Data Import -------------------------------------------------------------
rnb.set<-rnb.execute.import(data.source=list(data.dir, sample_sheet), data.type="idat.dir", verbose = TRUE)
#Check with Pavlo:s Error: In readChar(file, size, TRUE) : truncating string with embedded nuls

# analyze the bead counts
summarize.bead.counts<-function(bead.counts.M, bead.counts.U, method="min"){
  nrow.M<-nrow(bead.counts.M); ncol.M<-ncol(bead.counts.M)
  nrow.U<-nrow(bead.counts.U); ncol.U<-ncol(bead.counts.U)
  if(nrow.M!=nrow.U || ncol.M!=ncol.U){
    stop("Dimensions of bead count matrices differ")
  }else{
    bead.counts<-matrix(NA,nrow.M,ncol.M)
    rownames(bead.counts)<-rownames(bead.counts.M)
  }
  if(method=="min"){
    index1<-bead.counts.M<=bead.counts.U
    index1[is.na(index1)]<-FALSE
    bead.counts[index1]<-bead.counts.M[index1]
    index2<-bead.counts.U<bead.counts.M
    index2[is.na(index2)]<-FALSE
    bead.counts[index2]<-bead.counts.U[index2]
  }
  bead.counts
}

bead.matrix <- summarize.bead.counts(rnb.set@bead.counts.M,rnb.set@bead.counts.U)[,1:27]
apply(bead.matrix,2,median)%>%median
min(bead.matrix)
max(bead.matrix)
median(bead.matrix)
colSums(bead.matrix>3)%>%min/nrow(bead.matrix)*100


## Generate a density plot for the detection p-values
dpval <- data.frame(rnb.set@pval.sites)[,1:27]
colSums(dpval>0.001)%>%min

colnames(dpval) <- pheno(rnb.set)$Sample_Name[1:27]
dpval.melted <- reshape2::melt(dpval)
dpval.melted$ct <- stringr::str_split_fixed(dpval.melted$variable,pattern = "_",2)[,1]

colours <- c("lsk"="gray50",
             "cmp"="darkorange",
             "gmp"="darkgoldenrod2",
             "mep"="orangered4",
             "mono"="gold1",
             "neutro"="darksalmon",
             "bcell"="navy",
             "cd8"="lightsteelblue1",
             "cd4"="lightskyblue")

pdf(paste0(plot.dir,Sys.Date(),"_detection_pvals.pdf"))
ggplot(dpval.melted,aes(value,fill=ct))+
  #geom_histogram(bins=100)+
  theme_classic()+
  facet_wrap(~ct)+
  geom_bar(aes(y = (..count..)/784800)) + 
  scale_fill_manual(values=colours)+
  scale_y_continuous(labels = scales::percent)
ggplot(dpval.melted[which(dpval.melted$value>0.001),],aes(value,fill=ct))+
  #geom_histogram(bins=100)+
  theme_classic()+
  facet_wrap(~ct)+
  geom_bar(aes(y = (..count..)/784800)) + 
  scale_fill_manual(values=colours)+
  scale_y_continuous(labels = scales::percent)+
  ggtitle("Detection p-Values (>0.001) Percentage Total")
ggplot(dpval.melted[which(dpval.melted$value>0.001),],aes(value,fill=ct))+
  geom_histogram(bins=100)+
  theme_classic()+
  facet_wrap(~ct)+
  scale_fill_manual(values=colours)+
  ggtitle("Detection p-Values (>0.001) Count")+
  ylab("n (Sites)")+
  xlab("p-value")
dev.off()

## check if the sample annotation is correct
rnb.run.qc(rnb.set,dir.reports=paste0(report.dir,"qc_03_2021"))
rnb.set2 <- rnb.run.preprocessing(rnb.set,dir.reports=paste0(report.dir,Sys.Date(),"_preprocessing"))$rnb.set

#save the RnBeads Object
saveRDS(rnb.set2, paste0(data.dir,Sys.Date(),"_rnbeadsset_scaling.RDS"))
#save.rnb.set(rnb.set2, paste0(data.dir,Sys.Date(),"_rnbeadsset"),archive=T)

# save the lighter version of the RnBeads Object without p-values  --------
rnbeads.light <- as(rnb.set2, "RnBeadSet")
rnbeads.light@pval.sites <- NULL
rnbeads.light@covg.sites <- NULL
saveRDS(rnbeads.light, paste0(data.dir,Sys.Date(),"_rnbeadsset_betas_scaling.RDS"))
#save.rnb.set(rnbeads.light, paste0(data.dir,Sys.Date(),"_rnbeadsset_betas"),archive=T)
