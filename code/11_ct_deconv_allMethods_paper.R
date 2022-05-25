###################################################################
#
# Script: 11_ct_deconv_allMethods_paper.R
# Apply cell type deconvolution methods
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(EpiDISH)
library(MeDeCom)
library(Metrics)
library(ggpubr)
library(viridis)
library(gridExtra)
library(grid)
#remotes::install_github("lutsik/MeDeCom")

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
save.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/RDS/differential/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"

# Load the dataset --------------------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb.full@pheno$CT <- rnb.full@pheno$CellType
rnb.full@pheno$CT[28:36] <- NA
meth.full <- meth(rnb.full,row.names=T)
colnames(meth.full) <- rnb.full@pheno$Sample_Name

# Load the FACS data -----------------------------------------
facs <- read.delim(paste0(table.dir,"2021-05-19_cell_type_cont_red.txt"),stringsAsFactors = F,dec = ",")
rownames(facs) <- facs$Sample
facs <- facs[,-1]
colnames(facs) <- c("MyP","EryP","LSK","B-cell","Monocyte","Neutrophil","T-Cells")

# Meth CT -----------------------------------------------------------------
rnb.ct <- remove.samples(rnb.full,c(1:3,7:12,28:36))
rnb.ct@pheno$CT[1:3] <- "EryP"
rnb.ct@pheno$CT[4:9] <- "T-Cells"
rnb.ct@pheno$CT <- factor(rnb.ct@pheno$CT)

# Run Differential
cmp.cols <- "CT"
diffmeth <- rnb.execute.computeDiffMeth(rnb.ct,cmp.cols,region.types=c(),disk.dump = F) #regions turned off
comps <- get.comparisons(diffmeth)
levels(rnb.ct@pheno$CT)

## Select highly differential sites
n_site <- 50
diff.list <- list()
for(i in 1:length(comps)){
  comparison <- comps[i]
  tab.sites <- get.table(diffmeth, comparison, "sites", return.data.frame=TRUE)
  tab.sites$site <- rownames(meth.full)
  #ind1 <- which(tab.sites$min.g1>0.7&tab.sites$max.g2<0.2&tab.sites$diffmeth.p.adj.fdr<0.05&tab.sites$sd.g1<0.1&tab.sites$sd.g2<0.1)
  #message(paste0("Index1: ",length(ind1)))
  ind2 <- which(tab.sites$max.g1<0.2&tab.sites$min.g2>0.7&tab.sites$diffmeth.p.adj.fdr<0.05&tab.sites$sd.g1<0.1&tab.sites$sd.g2<0.1)
  message(paste0("Index2: ",length(ind2)))
  #ind <- c(ind1,ind2)
  ind <- ind2
  df <- tab.sites[ind,]
  if(nrow(df)>n_site){
    ind3 <- order(abs(df$mean.diff),decreasing = T)[1:n_site]
    df <- df[ind3,]
  }
  df$CellType <- rep(levels(rnb.ct@pheno$CT)[i],nrow(df))
  diff.list[[levels(rnb.ct@pheno$CT)[i]]] <- df
}

diff.df <- do.call(rbind,diff.list)
table(diff.df$CellType)
write.table(diff.df,paste0(table.dir,Sys.Date(),"_ct_deconv_cpgs.txt"),sep="\t",quote=F)

# Plot the differential sites ---------------------------------------------
filtered <- which(!rownames(meth.full)%in%diff.df$site)
rnb.sub <- remove.sites(rnb.full,filtered)
rnb.sub <- remove.samples(rnb.sub,c(1:3,7:12,32:36))
rnb.sub@pheno$CT[1:3] <- "EryP"
rnb.sub@pheno$CT[4:9] <- "T-Cells"
rnb.sub@pheno$CT <- factor(rnb.sub@pheno$CT)

plot.sub <- meth(rnb.sub,row.names=T)
colnames(plot.sub) <- pheno(rnb.sub)$Sample_Name

##colour scheme
annot.df <- rnb.sub@pheno[,c("Replicate","CT")]
rownames(annot.df) <- rnb.sub@pheno$Sample_Name

## define the colours
colours <- c("EryP"="orangered4",
             "Monocyte"="gold1",
             "Neutrophil"="darksalmon",
             "B-cell"="navy",
             "T-Cells"="lightskyblue")

anno.colour <- list("CT"=colours,
                    "Replicate"=c("gray90","gray60","gray20"))


pdf(paste0(plot.dir,Sys.Date(),"_rnbsub.pdf"))
pheatmap(plot.sub[,-c(19:22)][,c(4:9,16:18,13:15,1:3,10:12)],
         annotation_col = annot.df[-c(19:22),],
         annotation_colors = anno.colour,
         cluster_cols = F,
         show_rownames = F)
pheatmap(plot.sub,
         show_rownames = F)
dev.off()

# RnBeads Contributions --------------------------------------------------
ct.obj<-rnb.execute.ct.estimation(rnb.sub, cell.type.column="CT", test.max.markers=NULL, top.markers=sum(table(diff.df$CellType)))
plotting <- ct.obj$contributions[1:4,]
rownames(plotting) <- rnb.sub@pheno$Sample_Name[is.na(rnb.sub@pheno$CT)]
round(plotting,2)

# EpiDish -----------------------------------------------------------------

mat.full <- cbind("EryP"=apply(plot.sub[,1:3],1,mean),
                  "T-Cells"=apply(plot.sub[,4:9],1,mean),
                  "B-cell"=apply(plot.sub[,10:12],1,mean),
                  "Monocyte"=apply(plot.sub[,13:15],1,mean),
                  "Neutrophil"=apply(plot.sub[,16:18],1,mean))
preds.full <- plot.sub[,19:22]

# RPC Mode
ed.rpc <- epidish(beta.m = preds.full, ref.m = mat.full, method = "RPC") 
ed.rpcF <- ed.rpc$estF

# Cybersort Mode
ed.cbs <- epidish(beta.m = preds.full, ref.m = mat.full, method = "CBS") 
ed.cbsF <- ed.cbs$estF

# MeDeCom -----------------------------------------------------------------
#medecom.result<-runMeDeCom(plot.sub, Ks=2:10, lambdas=c(0,10^(-5:-1)),NINIT=10, NFOLDS=10, ITERMAX=100, NCORES=20)
#saveRDS(medecom.result,paste0(save.dir,Sys.Date(),"_medecom_allPops.RDS"))
medecom.result <- readRDS(paste0(save.dir,"2021-09-13_medecom_allPops.RDS"))

## check results
pdf(paste0(plot.dir,Sys.Date(),"_medecom_metric.pdf"))
plotParameters(medecom.result)
dev.off()

# best parameters are k=5 and lambda=0.01
pdf(paste0(plot.dir,Sys.Date(),"_medecom_lmcs.pdf"))
plotLMCs(medecom.result, K=5, lambda=0.01, type="dendrogram", Tref=mat.full, center=TRUE)
plotLMCs(medecom.result, K=5, lambda=0.01, type="heatmap", Tref=mat.full)
plotProportions(medecom.result, K=5, lambda=0.01, type="barplot")
dev.off()

prop<-getProportions(medecom.result, K=5, lambda=0.01)
colnames(prop) <- colnames(plot.sub)
rownames(prop) <- c("Neutrophil","Monocyte","T-Cells","B-cell","EryP")
medecom.stat <- t(prop[,19:22])

# Compare the metrics -----------------------------------------------------
facs <- facs[,c(2,4:7)]/100
epidish.stat.rpc <- ed.rpcF[rownames(facs),colnames(facs)]
epidish.stat.cbs <- ed.cbsF[rownames(facs),colnames(facs)]
rnbeads.stat <- plotting[rownames(facs),colnames(facs)]
medecom.stat <- medecom.stat[rownames(facs),colnames(facs)]

plotting.df <- rbind(reshape2::melt(as.matrix(epidish.stat.rpc)),
                     reshape2::melt(as.matrix(epidish.stat.cbs)),
                     reshape2::melt(as.matrix(rnbeads.stat)),
                     reshape2::melt(as.matrix(medecom.stat)))


plotting.df$tech <- rep(c("EpiDishRPC","EpiDishCBS","RnBeads","MeDeCom"),each=20)
facs.melt <- reshape2::melt(as.matrix(facs))
plotting.merged <- merge(plotting.df,facs.melt,by.x=c("Var1","Var2"),by.y=c("Var1","Var2"))

##sort the cells
plotting.merged$Var2 <- factor(plotting.merged$Var2,levels=c("B-cell","Monocyte","Neutrophil","T-Cells","EryP"))

## combined plot
pdf(paste0(plot.dir,Sys.Date(),"_ct_deconv_error.pdf"))
ggplot(plotting.merged,aes(value.x,value.y,shape=tech,colour=plotting.merged$Var2))+
  geom_abline(intercept = 0,alpha=.2,linetype = 2)+
  geom_point(size=2.5,alpha=.5)+
  theme_classic()+
  coord_cartesian(xlim = c(0, 0.5),ylim = c(0, 0.5))+
  scale_colour_manual(values=colours)+
  ylab("FACS [Fraction of CD45+]")+
  xlab("Predicted Fraction")
dev.off()

## individual plots
pdf(paste0(plot.dir,Sys.Date(),"_ct_deconv_ind.pdf"),width = 15,height = 5)
ggplot(plotting.merged,aes(value.x,value.y,colour=plotting.merged$Var2))+
  stat_cor(label.y.npc="top", label.x.npc = "left")+
  geom_abline(intercept = 0,alpha=.2,linetype = 2)+
  geom_point(aes(shape=Var1),size=4,alpha=.8)+
  facet_wrap(~tech,nrow=1)+
  theme_classic()+
  theme(legend.position="bottom")+
  coord_cartesian(xlim = c(0, 0.5),ylim = c(0, 0.5))+
  scale_colour_manual(values=colours)+
  ylab("FACS [Fraction of CD45+]")+ 
  xlab("Predicted Fraction")
dev.off()



# MAE ---------------------------------------------------------------------

list.stats <- list(facs,epidish.stat.rpc,epidish.stat.cbs,rnbeads.stat,medecom.stat)
list.mae <- lapply(list.stats,function(x){
for(i in 1:ncol(facs)){
  if(i==1){
    met <- mae(facs[,i],x[,i])
  }else{
    met <- c(met,mae(facs[,i],x[,i]))
  }
}
  return(met)})
mae.df <- do.call(rbind,list.mae)
rownames(mae.df) <- c("FACS","EpiDishRPC","EpiDishCBS","RnBeads","MeDeCom")
colnames(mae.df) <- colnames(facs)

pdf(paste0(plot.dir,Sys.Date(),"_mae_heatmap_all.pdf"))
pheatmap(mae.df[-1,],
         #breaks=seq(0, 0.3, by = 0.01),
         #color = colorRampPalette(c("steelblue", "white","darkred"))(length(seq(0, 0.3, by = 0.01))),
         color = colorRampPalette(c("steelblue", "white","darkred"))(100),
         #color=rev(viridis(100)),
         cluster_rows = F,
         cluster_cols = T,
         treeheight_col = 0,
         main = "MAE",
         display_numbers = T
        )
pheatmap(mae.df[-1,-1],
         #breaks=seq(0, 0.1, by = 0.01),
         color = colorRampPalette(c("steelblue", "white","darkred"))(100),
         #color = colorRampPalette(c("steelblue", "white","darkred"))(length(seq(0, 0.1, by = 0.01))),
         #color=rev(viridis(100)),
         cluster_rows = F,
         cluster_cols = T,
         treeheight_col = 0,
         main = "MAE",
         display_numbers = T
)
dev.off()


# JAK2 mice ---------------------------------------------------------------
plotting.df2 <- rbind(reshape2::melt(as.matrix(epidish.stat.rpc)),
                     reshape2::melt(as.matrix(epidish.stat.cbs)),
                     reshape2::melt(as.matrix(rnbeads.stat)),
                     reshape2::melt(as.matrix(medecom.stat)),
                     reshape2::melt(as.matrix(facs)))
plotting.df2$tech <- rep(c("EpiDishRPC","EpiDishCBS","RnBeads","MeDeCom","facs"),each=20)
plotting.df2$id <- stringr::str_split_fixed(plotting.df2$Var1,"_",2)[,1]
plotting.df2$id <- factor(plotting.df2$id,levels=c("jak2","cd45"))

pdf(paste0(plot.dir,Sys.Date(),"_jak2_cd45_cell_types.pdf"))
ggplot(plotting.df2[plotting.df2$Var2=="Monocyte",],aes(tech,value,fill=id,colour=id))+
  geom_bar(width = 0.8, position = position_dodge(width = 0.9), stat = "summary", fun = "mean")+
  geom_point(aes(shape=Var1),colour="black",size=3)+
  theme_classic()
dev.off()

