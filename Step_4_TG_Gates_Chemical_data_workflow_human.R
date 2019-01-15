#Author: Tatyana Doktorova
#Version:12_12_2018
# The script uses the information extracted from TG Gates and filtering of interesting genes according to fold change and p-value.
# Tge genes are re-mapped to Ensembl ID and NGTX-specific  genes extracted and visualized.


#BiocInstaller::biocLite('grimbough/biomaRt')
library(biomaRt)

library("biomaRt")
library("dplyr")
library(collapsibleTree)

#Generation of GTX list of genes

TG_Gates_human_log2fold_GTX<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/TGG_HCC_data_log2fold_human_GTX.txt", header=TRUE)
GTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/TGG_in_HCC_samples_GTX.txt", header=TRUE)
HCC_Toxcast_CTD<-read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/HCC_Toxcast_CTD.csv", header=TRUE)
Unique_genes = HCC_Toxcast_CTD %>% distinct(lhs)

#Filtering according to p value and fold change
GTX_filtering<-TG_Gates_human_log2fold_GTX[!is.na(TG_Gates_human_log2fold_GTX$pvalue),]
GTX_pval0.05<- GTX_filtering[GTX_filtering$pvalue<0.01,]
GTX_pval0.05_FC_up_1<-GTX_pval0.05[GTX_pval0.05$ value > 1,]
GTX_pval0.05_FC_down_1<-GTX_pval0.05[GTX_pval0.05$ value < -1,]
GTX_significant<- rbind(GTX_pval0.05_FC_down_1,GTX_pval0.05_FC_up_1)
GTX_significant_unique <- GTX_significant[!duplicated(GTX_significant[,2]),]


#Re-mapping to ensembl identifier
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

mapping <- getBM(attributes = c("affy_hg_u133_plus_2",
                                "ensembl_gene_id", "hgnc_symbol"), filters = "affy_hg_u133_plus_2",
                 values = GTX_significant_unique$assayId, mart = ensembl)

GTX_GL<-merge (mapping,GTX_significant_unique, by.x="affy_hg_u133_plus_2", by.y="assayId")

GTX_GL_final <- GTX_GL[!duplicated(GTX_GL[,3]),]
GTX_GL_final<-merge(GTX_GL_final,GTX_sampleID)

GTX_compounds_GL<-GTX_GL_final %>% distinct(TGG_compoundName)

#Generation of NGTX list of gene

TG_Gates_human_log2fold_NGTX<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/TGG_HCC_data_log2fold_human_NGTX.txt", header=TRUE)
NGTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/TGG_in_HCC_samples_NGTX.txt", header=TRUE)

NGTX_filtering<-TG_Gates_human_log2fold_NGTX[!is.na(TG_Gates_human_log2fold_NGTX$pvalue),]
NGTX_pval0.05<- NGTX_filtering[NGTX_filtering$pvalue<0.01,]
NGTX_pval0.05_FC_up_1<-NGTX_pval0.05[NGTX_pval0.05$ value > 1,]
NGTX_pval0.05_FC_down_1<-NGTX_pval0.05[NGTX_pval0.05$ value < -1,]
NGTX_significant<- rbind(NGTX_pval0.05_FC_down_1,NGTX_pval0.05_FC_up_1)


NGTX_significant_unique <- NGTX_significant[!duplicated(NGTX_significant[,2]),]



mapping_NGTX <- getBM(attributes = c("affy_hg_u133_plus_2",
                                "ensembl_gene_id", "hgnc_symbol"), filters = "affy_hg_u133_plus_2",
                 values = NGTX_significant_unique$assayId, mart = ensembl)

NGTX_GL<-merge (mapping_NGTX,NGTX_significant_unique, by.x="affy_hg_u133_plus_2", by.y="assayId")

NGTX_GL_final <- NGTX_GL[!duplicated(NGTX_GL[,3]),]
NGTX_GL_final<-merge(NGTX_GL_final,NGTX_sampleID)

NGTX_compounds_GL<-NGTX_GL_final %>% distinct(TGG_compoundName)


#Generation of NC list of gene

TG_Gates_human_log2fold_NC<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/TGG_HCC_data_log2fold_human_NC.txt", header=TRUE)
NC_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/TGG_in_HCC_samples_NC.txt", header=TRUE)

NC_filtering<-TG_Gates_human_log2fold_NC[!is.na(TG_Gates_human_log2fold_NC$pvalue),]
NC_pval0.05<- NC_filtering[NC_filtering$pvalue<0.01,]
NC_pval0.05_FC_up_1<-NC_pval0.05[NC_pval0.05$ value > 1,]
NC_pval0.05_FC_down_1<-NC_pval0.05[NC_pval0.05$ value < -1,]
NC_significant<- rbind(NC_pval0.05_FC_down_1,NC_pval0.05_FC_up_1)


NC_significant_unique <- NC_significant[!duplicated(NC_significant[,2]),]



mapping_NC <- getBM(attributes = c("affy_hg_u133_plus_2",
                                     "ensembl_gene_id", "hgnc_symbol"), filters = "affy_hg_u133_plus_2",
                      values = NC_significant_unique$assayId, mart = ensembl)

NC_GL<-merge (mapping_NC,NC_significant_unique, by.x="affy_hg_u133_plus_2", by.y="assayId")

NC_GL_final <- NC_GL[!duplicated(NC_GL[,3]),]
NC_GL_final<-merge(NC_GL_final,NC_sampleID)

NC_compounds_GL<-NC_GL_final %>% distinct(TGG_compoundName)

GTX_GL_genes<-unique(subset(GTX_GL_final, select=c("hgnc_symbol")))
NGTX_GL_genes<-unique(subset(NGTX_GL_final, select=c("hgnc_symbol")))
NC_GL_genes<-unique(subset(NC_GL_final, select=c("hgnc_symbol")))


Genes_common<-intersect(GTX_GL_genes, NGTX_GL_genes,NC_GL_genes)


NGTX_specific_1<-setdiff(NGTX_GL_genes, GTX_GL_genes)
NGTX_specific_2<-setdiff(NGTX_GL_genes, NC_GL_genes)
NGTX_specific<-merge(NGTX_specific_1, NGTX_specific_2)
colnames(Unique_genes) <- c("hgnc_symbol")
Unique_genes<-data.frame(Unique_genes)
NGTX_specific<-data.frame(NGTX_specific)

Genes_TGGATEs_CTD_Toxcast<-intersect(NGTX_specific, Unique_genes)

NGTX_induced_HCC<-merge(Genes_TGGATEs_CTD_Toxcast, HCC_Toxcast_CTD, by.x="hgnc_symbol", by.y="lhs")

collapsibleTree(df=NGTX_induced_HCC, c("Major", "Minor","Parent.Name", "hgnc_symbol"),  fill = "lightgreen")
collapsibleTree(df=HCC, c("Major", "Minor","Parent.Name", "lhs"),  fill = "lightgreen")

collapsibleTree(df=NGTX_induced_HCC, c("hgnc_symbol", "Major","Minor","Parent.Name"),  fill = "transparent")
