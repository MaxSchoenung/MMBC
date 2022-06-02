###################################################################
#
# Script: 10_dmp_methChanges_2021-05-23.R
# Methylation change in DMPs
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
library(dplyr)
library(pheatmap)
library(RnBeads)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"

# Load the HSC dmps -------------------------------------------------------
dmps.all <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
dmps <- dmps.all[dmps.all$mean.diff>0,]

# Load the methylation data -----------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)
rm(rnb.full)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name


# Changes Monocyte DMPs ---------------------------------------------------
list.change <- list()
for(i in names(table(dmps$CellType))){
ind <- dmps[dmps$CellType==i,"site"]
message(i)
beta.ind <- meth.rnb[ind,]
beta.melt <- reshape2::melt(beta.ind)
beta.melt$id <- stringr::str_split_fixed(beta.melt$Var2,"_",2)[,1]
beta.mean <- beta.melt%>%
  group_by(Var1,id)%>%
  summarise_at("value",mean)
beta.mean$ct <- rep(i,nrow(beta.mean))
list.change[[i]] <- beta.mean}
plot.df <- do.call(rbind,list.change)
head(plot.df)

## specify the colouring
colours <- c("lsk"="gray50",
             "cmp"="darkorange",
             "gmp"="darkgoldenrod2",
             "mep"="orangered4",
             "mono"="gold1",
             "neutro"="darksalmon",
             "bcell"="navy",
             "cd8"="lightsteelblue1",
             "cd4"="lightskyblue")


plot.df$id <- factor(plot.df$id,levels=c("lsk","bcell","cd8","cd4","mep","cmp","gmp","neutro","mono"))
plot.df$ct <- factor(plot.df$ct,levels=c("CD4_t-cell","CD8_t-cell","B-cell","MEP","CMP","GMP","Monocyte","Neutrophil"))
head(plot.df)

pdf(paste0(plot.dir,Sys.Date(),"_dmrs_mean_facet.pdf"),height =10)
ggplot(plot.df,aes(id,value,fill=id))+
  geom_boxplot()+
  facet_wrap(~ct,ncol=2)+
  theme_classic()+
  scale_fill_manual(values = colours)+
  ylab("Beta-Value")+
  xlab("Populations")+
  theme(legend.position="bottom")
dev.off()

# Delta Meth --------------------------------------------------------------
#first test with the mono dmps
dmps.mono <- dmps.all[dmps.all$CellType=="Monocyte",]
meth.mean <-   data.frame(
  "LSK"=rowMeans(meth.rnb[,1:3]),
  "CMP"=rowMeans(meth.rnb[,7:9]),
  "GMP"=rowMeans(meth.rnb[,10:12]),
  "Mono"=rowMeans(meth.rnb[,22:24]))
  #"Neutro"=rowMeans(meth.rnb[,25:27]))

meth.sub <- meth.mean[dmps.mono$site,]
meth.diff <- meth.sub-meth.sub$LSK
meth.diff <- cbind(meth.diff,"hypo"=as.factor(dmps.mono$mean.diff>0))
diff.melt <- reshape2::melt(meth.diff,id.vars="hypo")

pdf(paste0(plot.dir,Sys.Date(),"_dmrs_mono_diff_lsk.pdf"))
ggplot(diff.melt,aes(variable,value,fill=hypo))+
  geom_hline(yintercept=0)+
  geom_violin()+
  theme_classic()+
  scale_fill_manual(values = c("TRUE"="lightblue","FALSE"="firebrick"))
dev.off()


