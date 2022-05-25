###################################################################
#
# Script: 04_CT_tree.R
# Phylogenetic tree based on 5k mvCpGs
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
library(ape)

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

# 5k mv CpGs --------------------------------------------------------------
meth.sd <- apply(meth.rnb,1,sd)
sorted <- order(meth.sd,decreasing = T)[1:5000]
meth.var <- meth.rnb[sorted,]

# Calculate the distance metrix -------------------------------------------------------
dist.mat <- dist(t(meth.var),method = "manhattan")

# Minimal evolution -------------------------------------------------------
evol <- fastme.bal(dist.mat)
pdf(paste0(plot.dir,Sys.Date(),"_5k_mv_tree.pdf"))
plot.phylo(evol,
           "u")
dev.off()




