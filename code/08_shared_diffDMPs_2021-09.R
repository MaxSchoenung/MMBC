###################################################################
#
# Script: 08_clustering_diffDMPs.R
# Clustering of differentiation DMPs (diffDMPs)
# Date: 2021-05-14
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
#remotes::install_github("epigen/RnBeads",ref="feature/MouseMethylationBeadChip")
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Load the dataset --------------------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name
dim(meth.rnb)

# Load the HSC dmps -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)

# Hypo and hypermethylated DMPs per cell type which are unique and shared -----------------------------------------------------------
sites <- table(dmps$site) #how often does a site appear
attributes(sites)$class <- "matrix"
non_unique <- names(sites[sites>1]) #extract the names of the nun unique sites
dmps$non_unique <- as.numeric(dmps$site%in%non_unique) #add this to the DMP annotation
dmps$hypermethylated <- as.numeric(dmps$mean.diff<0) #which DMPs are hypermethylated
dmps.unique <- dmps[!duplicated(dmps$site),]
plot.df <- table(dmps$CellType,dmps$non_unique,dmps$hypermethylated) #hypermethylated DMPs per cell type
attributes(plot.df)$class <- "matrix"

# add annotation in hypo and hyper
plot.new <- cbind(plot.df[,,1],rep("hypo",nrow(plot.df[,,1])))
plot.new <- rbind(plot.new,cbind(plot.df[,,2],rep("hyper",nrow(plot.df[,,2]))))

#add the total number of DMPs
total.df <- table(dmps.unique$non_unique,dmps.unique$hypermethylated)
attributes(total.df)$class <- "matrix"
total.df <- cbind(t(total.df),c("hypo","hyper"))
rownames(total.df) <- rep("total",2)
plot.new <- rbind(plot.new,total.df)

## shape the data for plotting
colnames(plot.new) <- c("unique","shared","diff")
plot.new <- cbind(plot.new,"CellType"=rownames(plot.new))
plot.new <- as.data.frame(plot.new)
plot.melt <- reshape2::melt(plot.new,id.vars=c("diff","CellType"))
plot.melt$CellType <- factor(plot.melt$CellType,levels=rev(c("CMP","GMP","MEP","B-cell","Monocyte","Neutrophil","CD8_t-cell","CD4_t-cell","total")))
plot.melt$variable <- factor(plot.melt$variable)
plot.melt$value <- as.numeric(plot.melt$value)
plot.melt[plot.melt$diff=="hypo","value"] <- 0-plot.melt[plot.melt$diff=="hypo","value"]

# write to file
write.table(plot.melt,paste0(table.dir,Sys.Date(),"_dmps_hyper_hypo.txt"),sep="\t",row.names = T,quote = F,col.names = T)


## plot
pdf(paste0(plot.dir,Sys.Date(),"_hsc_dmps_number.pdf"))
ggplot(plot.melt,aes(value,CellType,fill=variable))+
  geom_bar(stat="identity")+
  theme_classic()+
  geom_vline(xintercept=0)+
  xlim(c(-32000,32000))
dev.off()

# Cell Type x Site --------------------------------------------------------
cxs <- table(dmps$site,dmps$CellType)
cxs.cor <- cor(cxs)

reshape_fun <- function(x){
  mat <- matrix(0,ncol=ncol(x),nrow=ncol(x))
  colnames(mat) <- colnames(x)
  rownames(mat) <- colnames(x)
  for(i in 1:nrow(x)){
    if(sum(x[i,])>1){
    ind <- which(x[i,]==1)
    mat[ind,ind] <- mat[ind,ind]+1}
  }
  return(mat)
}


shared <- reshape_fun(cxs)
sort.vec <- c("CMP","GMP","Monocyte","Neutrophil","MEP","B-cell","CD8_t-cell","CD4_t-cell")

write.table(shared[sort.vec,sort.vec],paste0(table.dir,Sys.Date(),"shared_DMPs.txt"),sep="\t",row.names = T,quote = F,col.names = T)

shared <- apply(shared,2,function(x){(x/c(table(dmps$CellType[dmps$non_unique==1])))*100}) #normalize by total DMPs
write.table(shared[sort.vec,sort.vec],paste0(table.dir,Sys.Date(),"shared_DMPs_perc.txt"),sep="\t",row.names = T,quote = F,col.names = T)

pdf(paste0(plot.dir,Sys.Date(),"_CellType_x_Site.pdf"))
pheatmap(shared[sort.vec,sort.vec],
         show_rownames = T,
         treeheight_col = 0,
         treeheight_row = 0,
         cluster_cols = F,
         cluster_rows = F,
         na_col="white",
         color = colorRampPalette(c("white","darkred"))(100),
         main = "Percentage Shared per LSK-CellType DMPs (Row-wise)"
         )
dev.off()
