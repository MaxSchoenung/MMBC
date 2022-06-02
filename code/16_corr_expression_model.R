###################################################################
#
# Script: 16_corr_expression_model.R
# Correlation DMPs and expression
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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DESeq2)
library(ChIPpeakAnno)
library(ggplot2)
library(VennDiagram)
library(caret)
library(ggpointdensity)
library(viridis)
library(foreach)
library(parallel)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"

# Load the methylation data -----------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)
rm(rnb.full)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name

# Load the HSC dmps and novel CRE programs -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
programs <- read.delim(paste0(table.dir,"2021-07-30_new_cre_9_cluster.txt"),stringsAsFactors = F)

# Load the annotation -----------------------------------------------------
annot.full <- read.delim(paste0(table.dir,"2022-05-25_MMBC_annotation_modified.bed"),stringsAsFactors = F)
rownames(annot.full) <- annot.full$IlmnID
annot.full <- annot.full[annot.full$Chromosome%in%paste0("chr",1:19),]
annot <- annot.full[,1:3]
annot.ranges <- GenomicRanges::makeGRangesFromDataFrame(annot)
seqlevelsStyle(annot.ranges) = "UCSC"

## annotate the dmps
dmps.annot <- cbind(annot[dmps$site,],dmps)
programs.annot <- cbind(annot[rownames(programs),],programs)
dmp.ranges <- makeGRangesFromDataFrame(dmps.annot,keep.extra.columns = T)
programs.ranges <- makeGRangesFromDataFrame(programs.annot,keep.extra.columns = T)

# Load the RNA-Seq data ---------------------------------------------------
dds <- readRDS(paste0(data.dir,"haemopedia_dds_mmbc.RDS"))
dds <- DESeq(dds)
rld <- vst(dds,blind = F)
mat <- assay(rld)

# Annotate putative methCREs to nearest genes ----------------------------------------------------------------
mm10 <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="http://nov2020.archive.ensembl.org") #biomart for annotation
tx.id <- getBM(attributes=c("ensembl_gene_id","external_gene_name","ensembl_transcript_id_version","chromosome_name","transcription_start_site"),filters="ensembl_gene_id",values=rownames(counts(dds)),mart=mm10)
tx.id <- tx.id[tx.id$chromosome_name%in%1:19,]
tx.ranges <- makeGRangesFromDataFrame(tx.id, seqnames.field = "chromosome_name",start.field = "transcription_start_site",end.field = "transcription_start_site",keep.extra.columns = T,ignore.strand = T)
seqlevelsStyle(tx.ranges) <- "UCSC"

# Mean Methylation per grop -----------------------------------------------
mean.meth <- data.frame(
  "LSK"=rowMeans(meth.rnb[,1:3]),
  "MEP"=rowMeans(meth.rnb[,4:6]),
  "CMP"=rowMeans(meth.rnb[,7:9]),
  "GMP"=rowMeans(meth.rnb[,10:12]),
  "CD4"=rowMeans(meth.rnb[,13:15]),
  "CD8"=rowMeans(meth.rnb[,16:18]),
  "BCell"=rowMeans(meth.rnb[,19:21]),
  "Mono"=rowMeans(meth.rnb[,22:24]),
  "Neutro"=rowMeans(meth.rnb[,25:27])
)
mean.exp <- data.frame(
  "LSK"=rowMeans(mat[,15:16]),
  "MEP"=rowMeans(mat[,9:10]),
  "CMP"=rowMeans(mat[,23:24]),
  "GMP"=rowMeans(mat[,25:26]),
  "CD4"=rowMeans(mat[,27:28]),
  "CD8"=rowMeans(mat[,29:30]),
  "BCell"=rowMeans(mat[,1:8]),
  "Mono"=rowMeans(mat[,11:14]),
  "Neutro"=rowMeans(mat[,17:22])
)


# For each CpG in the subset select highly correlated genes ---------------
dist.var=1000000 #1Mb
dmp.ext <- programs.ranges+dist.var
overl.dmps <- findOverlaps(dmp.ext,tx.ranges)
overl.dmps.df <- as.data.frame(overl.dmps)

##filter to remove double associations and speed up
overl.dmps.df$id <- paste0(names(dmp.ext)[overl.dmps@from],"_",tx.ranges$ensembl_gene_id[overl.dmps@to])
filter.ind <- which(duplicated(overl.dmps.df$id))
overl.dmps.filt <- overl.dmps[-filter.ind]

corr_fun <- function(x){
  meth.vec <- as.numeric(mean.meth[names(dmp.ext)[x@from],])
  exp.vec <- as.numeric(mean.exp[tx.ranges$ensembl_gene_id[x@to],])
  cor.sub <- cor.test(meth.vec,exp.vec)
  mod.sub <- lm(exp.vec~meth.vec)
  data.frame("cpg"=names(dmp.ext)[x@from],
             "gene"=tx.ranges$ensembl_gene_id[x@to],
             "gene_id"=tx.ranges$external_gene_name[x@to],
             "cor_value"=cor.sub$estimate,
             "p_value"=cor.sub$p.value,
             "slope"=mod.sub$coefficients["meth.vec"],
             "dist_TSS"=tx.ranges[x@to]@ranges@start-programs.ranges[x@from]@ranges@start,
             "chr"=as.character(seqnames(tx.ranges[x@to])),
             "TSS_start"=tx.ranges[x@to]@ranges@start,
             "cpgs_start"=programs.ranges[x@from]@ranges@start
             )
}


message(paste0(Sys.time()," - Start Splitting dataset with ",length(overl.dmps.filt)," Matches"))
overl.list = lapply(as.list(1:length(overl.dmps.filt)), function(x) overl.dmps.filt[x,])
message(paste0(Sys.time()," - Start Correlation Calculation"))
corr.list  <- mclapply(overl.list,corr_fun,mc.cores=10)
message(paste0(Sys.time()," - Start Binding the List"))
cor.cpg.gene <- do.call(rbind,corr.list)
message(paste0(Sys.time()," - Done!!"))
nrow(cor.cpg.gene)==length(overl.dmps.filt)

# P-Value correction and save ---------------------------------------
work.cor <- cor.cpg.gene
work.cor$ind <- paste0(cor.cpg.gene$cpg,"_",cor.cpg.gene$gene)
cor.distinct <- work.cor[!duplicated(work.cor$ind),]
cor.distinct$p_adj <- p.adjust(cor.distinct$p_value,method = "BH")
max(table(cor.distinct$cpg))

## analyze the distribution of significant correlating CpGs
cor.sig <- cor.distinct[cor.distinct$p_adj<0.01&abs(cor.distinct$slope)>4,]

saveRDS(cor.distinct,paste0(data.dir,Sys.Date(),"_novelCRE_9programs_cpg_gene_corr.RDS"))
#cor.distinct <- readRDS(paste0(data.dir,"2021-09-17_novelCRE_9programs_cpg_gene_corr.RDS"))
write.table(cor.distinct,paste0(table.dir,Sys.Date(),"_novelCRE_9programs_cpg_gene_corr.txt"),sep="\t",quote = F,row.names = F)
saveRDS(cor.sig,paste0(data.dir,Sys.Date(),"_novelCRE_9programs_cpg_gene_corr_sig_only.RDS"))
write.table(cor.sig,paste0(table.dir,Sys.Date(),"_novelCRE_9programs_cpg_gene_corr_sig_only.txt"),sep="\t",quote = F,row.names = F)

test.file <- read.delim(paste0(table.dir,"2021-08-02_novelCRE_9programs_cpg_gene_corr.txt"))
all(round(test.file$p_value,8)==round(cor.distinct$p_value,8))

dist.vec <- table(cor.sig$cpg)
dist.df <- data.frame(
  "number"=c(sum(dist.vec==1),sum(dist.vec>=2&dist.vec<=5),sum(dist.vec>=6&dist.vec<=10),sum(dist.vec>=11&dist.vec<=20),sum(dist.vec>20)),
  "what"=factor(c("1","2-5","6-10","11-20",">20"),levels=c("1","2-5","6-10","11-20",">20")))
write.table(dist.df,paste0(table.dir,Sys.Date(),"_n_CRE_gene_assoc.txt"),sep="\t",quote = F,row.names = F)

pdf(paste0(plot.dir,Sys.Date(),"_novelCREs_TSS_distribution_sig_corrs.pdf"))
ggplot(dist.df,aes(what,number))+
  geom_bar(stat="identity")+
  theme_classic()+
  xlab("Number of associated genes")+
  ylab("#CpGs")
ggplot(cor.sig,aes(dist_TSS/1000))+
  geom_histogram(bins = 100)+
  theme_classic()+
  xlab("Distance to TSS in kb (sig. CpGs within 1Mb)")
dev.off()

pdf(paste0(plot.dir,Sys.Date(),"novelCREs_cor_distribution.pdf"))
ggplot(cor.distinct,aes(cor_value))+
  geom_density()+
  theme_classic()+
  xlab("Correlation Value (all DMPs within 1Mb to TSS)")
ggplot(cor.distinct[abs(cor.distinct$dist_TSS)<=1000,],aes(cor_value))+
  geom_density()+
  theme_classic()+
  xlab("Correlation Value (all DMPs within 1kb to TSS)")
ggplot(cor.distinct,aes(cor_value,colour=p_adj<0.05))+
  geom_density()+
  theme_classic()+
  xlab("Correlation Value (all DMPs within 1Mb to TSS)")
ggplot(cor.distinct,aes(cor_value,colour=p_adj<0.01))+
  geom_density()+
  theme_classic()+
  xlab("Correlation Value (all DMPs within 1Mb to TSS)")
dev.off()


# Plot the Thy1 gene as an example
thy1.plot <- cor.distinct[cor.distinct$gene_id=="Thy1",]
pdf(paste0(plot.dir,Sys.Date(),"_Thy1_contacts.pdf"),height = 4)
ggplot(thy1.plot,aes(dist_TSS,-log10(p_adj),fill=p_adj<0.01,colour=p_adj<0.01))+
  geom_bar(stat="identity",width=15)+
  theme_classic()+
  xlim(-1000000,1000000)
ggplot(thy1.plot,aes(dist_TSS,-log10(p_adj),fill=p_adj<0.01&abs(slope)>4,colour=p_adj<0.01&abs(slope)>4))+
  geom_bar(stat="identity",width=15)+
  theme_classic()+
  xlim(-1000000,1000000)
dev.off()

thy1.plot[thy1.plot$p_adj<0.01&abs(thy1.plot$slope)>4,]%>%dim
