###################################################################
#
# Script: 07_differential_LSK.R
# Differential calling from LSK to differentiated cell types
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

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"


# Load the dataset --------------------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)
rm(rnb.full)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name
dim(meth.rnb)

# Remove PTPN11 and WBs ---------------------------------------------------
rnb@pheno$CellType <- factor(rnb@pheno$CellType,
                             levels = c("LSK","CMP","GMP","MEP","Monocyte","Neutrophil","B-cell","CD8_t-cell","CD4_t-cell"))

# Run Differential  --------------------------------------------------------
cmp.cols <- "CellType"
diffmeth <- rnb.execute.computeDiffMeth(rnb,pheno.cols=cmp.cols,pheno.cols.all.pairwise="CellType",region.types=c(),disk.dump = F) #regions turned off
comps <- get.comparisons(diffmeth)[which(stringr::str_detect(get.comparisons(diffmeth),"LSK"))]

## based on p<0.05 and delta meth > 0.2
diff.list <- list()
for(i in 1:length(comps)){
comparison <- comps[i]
tab.sites <- get.table(diffmeth, comparison, "sites", return.data.frame=TRUE)
tab.sites$site <- rownames(meth.rnb)
ind <- which(abs(tab.sites$mean.diff)>0.2&tab.sites$diffmeth.p.adj.fdr<0.05)
df <- tab.sites[ind,]
df$CellType <- rep(levels(rnb@pheno$CellType)[i+1],nrow(df))
diff.list[[levels(rnb@pheno$CellType)[i+1]]] <- df
}

## bind the list
diff.df <- do.call(rbind,diff.list)
table(diff.df$CellType)
dim(diff.df)
unique(diff.df$site)%>%length()

## save the list
write.table(diff.df,paste0(table.dir,Sys.Date(),"_diff_meth_lsk_p_0.05_delt_0.2.txt"),sep="\t",quote=F)

##check with old list
old <- read.delim(paste0(table.dir,"2021-03-16_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
all(old$site==diff.df$site) #check!
dim(old)
dim(diff.df)
