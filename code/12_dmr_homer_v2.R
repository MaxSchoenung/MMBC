###################################################################
#
# Script: 12_dmr_homer_v2.R
# Run homer to analyze TF binding sites within DMPs
# Date: 2021-03-15
# Author: Maximilian Schönung
#
###################################################################

# Load Libraries ----------------------------------------------------------
library(dplyr)
library(pheatmap)
library(RnBeads)
library(RnBeads.mm10)
library(ggplot2)

# Set Directories ---------------------------------------------------------
setwd("/icgc/analysis/OE0565_projects/mmbc_hierarchy/")
data.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/"
plot.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/plots/"
table.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/"
homer.dir <- "/icgc/analysis/OE0565_projects/mmbc_hierarchy/tables/homer/"

# Load the HSC dmps -------------------------------------------------------
dmps <- read.delim(paste0(table.dir,"2021-05-14_diff_meth_lsk_p_0.05_delt_0.2.txt"),stringsAsFactors = F)
dmps <- dmps[dmps$mean.diff>0,] # focus on hypomethylated

# Load the methylation data -----------------------------------------------
rnb.full <- readRDS("/icgc/analysis/OE0565_projects/mmbc_hierarchy/data/2021-03-15_rnbeadsset_betas_scaling.RDS")
rnb <- remove.samples(rnb.full,28:36)
rm(rnb.full)

## access the meth and pheno data
pheno.rnb <- pheno(rnb)
pheno.rnb$Sample_Name
meth.rnb <- meth(rnb,row.names=T)
colnames(meth.rnb) <- pheno.rnb$Sample_Name

# Load MMBC anno ----------------------------------------------------------
sites=rnb.get.annotation("probesMMBC","mm10")
array.probes <- unlist(sites)
array.df <- as.data.frame(array.probes)
rownames(array.df) <- array.df$ID
array.df <- array.df[rownames(array.df)%in%rownames(meth.rnb),]

## shape for Homer
homer.df <- data.frame(array.df$seqnames,
                       array.df$start,
                       array.df$end,
                       array.df$ID,
                       rep("",nrow(array.df)),
                       array.df$strand,
                      stringsAsFactors = F)
#write.table(homer.df,paste0(homer.dir,Sys.Date(),"_homer_bg.bed"),row.names = F,col.names = F,sep="\t",quote = F)

# Export the CT DMPs for Homer --------------------------------------------
result.list <- list()
for(i in names(table(dmps$CellType))){
  message(paste0("Processing: ",i))
  system(paste0("mkdir ",homer.dir,i))
  input.data <- paste0(homer.dir,i,"/",Sys.Date(),"_homer_dmps_",i,".bed")
  bg.data <- paste0(homer.dir,i,"/",Sys.Date(),"_homer_bg_",i,".bed")
  write.table(homer.df[homer.df$array.df.ID%in%dmps[dmps$CellType==i,"site"],],input.data,row.names = F,col.names = F,sep="\t",quote = F)
  write.table(homer.df[!homer.df$array.df.ID%in%dmps[dmps$CellType==i,"site"],],bg.data,row.names = F,col.names = F,sep="\t",quote = F)
  system(paste0("/ngs_share/tools/homer/bin/findMotifsGenome.pl ",input.data,
               " mm10 ",homer.dir,i," -size -50,50 -bg ",bg.data," -nomotif"))
  tab <- read.delim(paste0(homer.dir,i,"/knownResults.txt"))
  colnames(tab) = c("motif_name", "consensus", "p_value", "log_p_value", "q_value", "n_targets", "perc_targets", "n_bg_targets", "perc_bg_targets")
  tab$CellType <- rep(i,nrow(tab))
  result.list[[i]] <- tab
}


# To Read -----------------------------------------------------------------
# result.list <- list()
# for(i in names(table(dmps$CellType))){
#   if(file.exists(paste0(homer.dir,i,"/knownResults.txt"))){
#     message("Found the respective homer results")
#   }
#   tab <- read.delim(paste0(homer.dir,i,"/knownResults.txt"))
#   colnames(tab) = c("motif_name", "consensus", "p_value", "log_p_value", "q_value", "n_targets", "perc_targets", "n_bg_targets", "perc_bg_targets")
#   tab$CellType <- rep(i,nrow(tab))
#   result.list[[i]] <- tab
# }

## concatenate table and write to file
result.df <- do.call(rbind,result.list)
write.table(result.df,paste0(table.dir,Sys.Date(),"_homer_results_raw_concatenated.txt"),sep="\t",row.names = F,quote = F,col.names = T)


# Analysis ----------------------------------------------------------------
result.df$enrichment_fc <- as.numeric(stringr::str_remove_all(result.df$perc_targets,"%"))/as.numeric(stringr::str_remove_all(result.df$perc_bg_targets,"%"))
result.df$tf <- toupper(stringr::str_split_fixed(result.df$motif_name,"\\(",2)[,1])
result.df$sig <- as.numeric(result.df$p_value<0.00001)

selected2 <- c("PAX5","EBF","EBF1","EBF2","TCFL2","TCF3","TCF7","LEF1","SLUG",
              "BATF","AP-1","JUNB","JUN-AP1","FOSL2",
              "ATF3","IRF1","IRF3","IRF8",
              "ATF1","ATF4","CHOP","EHF","ETV4","ETV1","ETV6",
              "NFIL3","PU.1","CEBP",
              "GATA6","GATA3","GATA2","GATA1",
              "GATA:SCL","KLF6","KLF5","KLF4","KLF3")


result.plot <- result.df[result.df$tf%in%selected2,c("CellType","p_value","tf","enrichment_fc")]
result.plot$CellType <- factor(result.plot$CellType,levels = c("MEP","CMP","GMP","Neutrophil","Monocyte","CD8_t-cell","CD4_t-cell","B-cell"))
key <- paste0(result.plot$CellType,"_",result.plot$tf)
result.plot <- result.plot[!duplicated(key),]
result.plot$neg_log_p <- -log10(result.plot$p_value)
result.plot$tf <- factor(result.plot$tf,levels=selected2)


pdf(paste0(plot.dir,Sys.Date(),"_homer_point_plot.pdf"),height = 8)
ggplot(result.plot,aes(CellType,tf,size=enrichment_fc,fill=neg_log_p>15))+
  geom_point(shape=21,colour="black",stroke=.5)+
  scale_fill_manual(values=c("white", "black"),name="Significant:\n-log10(p)>15")+
  theme_classic()+
  theme(line = element_blank())
dev.off()

