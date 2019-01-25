####################################################################
####Author: Noffisat Oki, Tatyana Doktorova
####Version: 1.0 08-28-2018
####Description: GTX Subset script:
####             The script uses the information extracted from TG Gates and 
###              filtering of interesting genes according to fold change and p-value.
####             Subsets chemical-gene data to filter for only genes active
####             in at least 3 chemicals 
###              Reannotation to Ensemb Identifier and visualizations
####             
####Notes: The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 20
####       should be used for this purpose.
####Potential issues:

##The root location folder should be specified here. 
#This is the main folder where any subfolders from which data will be read from or written to are located
Root_location <- "C:\\Users\\tatyana\\Documents\\Projects\\Advance project\\Advance"

##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

# Install libraries
library(data.table)
library(reshape2)
library("biomaRt")
library("dplyr")
library(collapsibleTree)

#Generation of GTX list of genes


#Generation of GTX list of genes for rat in vivo
#The files can be derived either from the previous Step_3B1, Step_3B2, Step3_B3 or directly downloaded from this link: https://mega.nz/#F!zDBmiKyZ

TG_Gates_rat_log2fold_GTX<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/rat HCC/TGG_data_log2fold_HCC_rat_GTX.txt", header=TRUE)
GTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Projects/Advance project/rat HCC/TGG_in_HCC_samples_rat_GTX.txt", header=TRUE)
HCC_Toxcast_CTD<-read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/HCC_Toxcast_CTD.csv", header=TRUE)
Unique_genes = HCC_Toxcast_CTD %>% distinct(lhs)

GTX_all<-merge(TG_Gates_rat_log2fold_GTX,GTX_sampleID, by.x="sampleId",by.y="sampleId")

#Filtering according to p value and fold change


GTX_filtering<-GTX_all[!is.na(GTX_all$pvalue),]

GTX_pval0.05<- GTX_filtering[GTX_filtering$pvalue<0.05,]
GTX_pval0.05_FC_up_1<-GTX_pval0.05[GTX_pval0.05$ value > 0.58,]
GTX_pval0.05_FC_down_1<-GTX_pval0.05[GTX_pval0.05$ value < -0.58,]
GTX_specific<- rbind(GTX_pval0.05_FC_down_1,GTX_pval0.05_FC_up_1)
#write.csv(GTX_significant, file="GTX_specific_rat.csv")

#Generation of NGTX list of genes
TG_Gates_rat_log2fold_NGTX<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_HCC_data_log2fold_rat_NGTX.txt", header=TRUE)
NGTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_in_HCC_samples_rat_NGTX.txt", header=TRUE)

NGTX_all<-merge(TG_Gates_rat_log2fold_NGTX,NGTX_sampleID, by.x="sampleId",by.y="sampleId")

NGTX_filtering<-NGTX_all[!is.na(NGTX_all$pvalue),]

NGTX_pval0.05<- NGTX_filtering[NGTX_filtering$pvalue<0.05,]
NGTX_pval0.05_FC_up_1<-NGTX_pval0.05[NGTX_pval0.05$ value > 0.58,]
NGTX_pval0.05_FC_down_1<-NGTX_pval0.05[NGTX_pval0.05$ value < -0.58,]
NGTX_specific<- rbind(NGTX_pval0.05_FC_down_1,NGTX_pval0.05_FC_up_1)
#write.csv(NGTX_significant, file="NGTX_specific.csv")

#Generation of NC list of genes
TG_Gates_rat_log2fold_NC<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_HCC_data_log2fold_rat_NC.txt", header=TRUE)
NC_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_in_HCC_samples_rat_NC.txt", header=TRUE)
#merge the two sources and remove NA values
NC_all<-merge(TG_Gates_rat_log2fold_NC,NC_sampleID, by.x="sampleId",by.y="sampleId")

NC_filtering<-NC_all[!is.na(NC_all$pvalue),]

NC_pval0.05<- NC_filtering[NC_filtering$pvalue<0.05,]
NC_pval0.05_FC_up_1<-NC_pval0.05[NC_pval0.05$ value > 0.58,]
NC_pval0.05_FC_down_1<-NC_pval0.05[NC_pval0.05$ value < -0.58,]
NC_specific<- rbind(NC_pval0.05_FC_down_1,NC_pval0.05_FC_up_1)
#write.csv(NC_significant, file="NC_specific.csv")



################################# GTX _2/3 chemicals with same gene data ##################
GTX_chemgene_data.long <- as.data.frame (GTX_specific)
##reformatting the data
GTX_chemgene_data.wide <- dcast(GTX_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_GTX_chemgene_data.wide <- (GTX_chemgene_data.wide)
recoding_GTX_chemgene_data.wide[,2:ncol(recoding_GTX_chemgene_data.wide)][recoding_GTX_chemgene_data.wide[,2:ncol(recoding_GTX_chemgene_data.wide)]>1] <- 1
rownames(recoding_GTX_chemgene_data.wide) <- recoding_GTX_chemgene_data.wide[,1]
recoding_GTX_chemgene_data.wide <- ((recoding_GTX_chemgene_data.wide[,-1]))
#write.table(recoding_GTX_chemgene_data.wide,file=file.path(outDir, "GTX_chem_gene_data_wide.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

Active_GTX_dataset<- recoding_GTX_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_GTX_dataset<- Active_GTX_dataset[,!apply((Active_GTX_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 2)]
chemicals <- as.data.frame(rownames(Active_GTX_dataset))
colnames(chemicals) <- "chemical"
Probes <- as.data.frame(colnames(Active_GTX_dataset))
Active_GTX_dataset <- cbind(chemicals, Active_GTX_dataset)
Active_GTX_dataset_T<-t(Active_GTX_dataset)
Active_GTX_dataset<-cbind(Row.Names=rownames(Active_GTX_dataset_T),Active_GTX_dataset_T)

write.table(Active_GTX_dataset,file=file.path(outDir, "Active_GTX_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(Probes,file=file.path(outDir, "Active_GTX_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)


ensembl <- useMart("ensembl")
ensembl = useDataset("rnorvegicus_gene_ensembl",mart=ensembl)

mapping_GTX <- getBM(attributes = c("affy_rat230_2",
                                    "ensembl_gene_id", "rgd_symbol"), filters = "affy_rat230_2",
                     values = Probes, mart = ensembl)


GTX_GL_final <- mapping_GTX[!duplicated(mapping_GTX[,3]),] 


################################# GTX _2/3 chemicals with same gene data##################
NGTX_chemgene_data.long <- as.data.frame (NGTX_specific)
##reformatting the data
NGTX_chemgene_data.wide <- dcast(NGTX_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_NGTX_chemgene_data.wide <- (NGTX_chemgene_data.wide)
recoding_NGTX_chemgene_data.wide[,2:ncol(recoding_NGTX_chemgene_data.wide)][recoding_NGTX_chemgene_data.wide[,2:ncol(recoding_NGTX_chemgene_data.wide)]>1] <- 1
rownames(recoding_NGTX_chemgene_data.wide) <- recoding_NGTX_chemgene_data.wide[,1]
recoding_NGTX_chemgene_data.wide <- ((recoding_NGTX_chemgene_data.wide[,-1]))
#write.table(recoding_GTX_chemgene_data.wide,file=file.path(outDir, "GTX_chem_gene_data_wide.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

Active_NGTX_dataset<- recoding_NGTX_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_NGTX_dataset<- Active_NGTX_dataset[,!apply((Active_NGTX_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 2)]
chemicals <- as.data.frame(rownames(Active_NGTX_dataset))
colnames(chemicals) <- "chemical"
NGTX_Probes <- as.data.frame(colnames(Active_NGTX_dataset))
Active_NGTX_dataset <- cbind(chemicals, Active_NGTX_dataset)
Active_NGTX_dataset_T<-t(Active_NGTX_dataset)
Active_NGTX_dataset<-cbind(Row.Names=rownames(Active_NGTX_dataset_T),Active_NGTX_dataset_T)

write.table(Active_NGTX_dataset,file=file.path(outDir, "Active_NGTX_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(NGTX_Probes,file=file.path(outDir, "Active_NGTX_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)




mapping_NGTX <- getBM(attributes = c("affy_rat230_2",
                                     "ensembl_gene_id", "rgd_symbol"), filters = "affy_rat230_2",
                      values = NGTX_Probes, mart = ensembl)


NGTX_GL_final <- mapping_NGTX[!duplicated(mapping_NGTX[,3]),] 



################################# NC Chemical X Gene data ##################
NC_chemgene_data.long <- as.data.frame (NC_specific)
##reformatting the data
NC_chemgene_data.wide <- dcast(NC_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_NC_chemgene_data.wide <- (NC_chemgene_data.wide)
recoding_NC_chemgene_data.wide[,2:ncol(recoding_NC_chemgene_data.wide)][recoding_NC_chemgene_data.wide[,2:ncol(recoding_NC_chemgene_data.wide)]>1] <- 1
rownames(recoding_NC_chemgene_data.wide) <- recoding_NC_chemgene_data.wide[,1]
recoding_NC_chemgene_data.wide <- ((recoding_NC_chemgene_data.wide[,-1]))
#write.table(recoding_GTX_chemgene_data.wide,file=file.path(outDir, "GTX_chem_gene_data_wide.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

Active_NC_dataset<- recoding_NC_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_NC_dataset<- Active_NC_dataset[,!apply((Active_NC_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 2)]
chemicals <- as.data.frame(rownames(Active_NC_dataset))
colnames(chemicals) <- "chemical"
NC_Probes <- as.data.frame(colnames(Active_NC_dataset))
Active_NC_dataset <- cbind(chemicals, Active_NC_dataset)
Active_NC_dataset_T<-t(Active_NC_dataset)
Active_NC_dataset<-cbind(Row.Names=rownames(Active_NC_dataset_T),Active_NC_dataset_T)

write.table(Active_NC_dataset,file=file.path(outDir, "Active_NC_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(NC_Probes,file=file.path(outDir, "Active_NC_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)




mapping_NC <- getBM(attributes = c("affy_rat230_2",
                                   "ensembl_gene_id", "rgd_symbol"), filters = "affy_rat230_2",
                    values = NC_Probes, mart = ensembl)


NC_GL_final <- mapping_NC[!duplicated(mapping_NC[,3]),] 

######Unique gene list indentification
GTX_GL_genes<-unique(subset(GTX_GL_final, select=c("rgd_symbol")))
NGTX_GL_genes<-unique(subset(NGTX_GL_final, select=c("rgd_symbol")))
NC_GL_genes<-unique(subset(NC_GL_final, select=c("rgd_symbol")))

##Common genees between GTX, NGTX, NC
Genes_common<-intersect(GTX_GL_genes, NGTX_GL_genes,NC_GL_genes)

## NGTX specific genes
NGTX_specific_1<-setdiff(NGTX_GL_genes, GTX_GL_genes)
NGTX_specific_2<-setdiff(NGTX_GL_genes, NC_GL_genes)
NGTX_specific<-merge(NGTX_specific_1, NGTX_specific_2)

NGTX_specific<-data.frame(NGTX_specific)

##Overlap between the ToxCast_CTD HCC-specific genes and TG Gates

Genes_TGGATEs_CTD_Toxcast<-intersect(NGTX_specific, Unique_genes)


NGTX_induced_HCC_overlap<-merge(Genes_TGGATEs_CTD_Toxcast, HCC_Toxcast_CTD, by.x="rgd_symbol", by.y="lhs")

NGTX<-rbind(NGTX_induced_HCC,NGTX_induced_HCC_overlap)

collapsibleTree(df=NGTX, c("rgd_symbol", "Major","Minor","Parent.Name"),  fill = "transparent")

collapsibleTree(df=NGTX_induced_HCC_overlap, c("Major", "Minor","Parent.Name", "rdg_symbol"),  fill = "lightgreen")


collapsibleTree(df=NGTX_induced_HCC_overlap, c("rgd_symbol", "Parent.Name"),  fill = "transparent")

collapsibleTree(df=NGTX_induced_HCC_overlap, c("rdg_symbol", "Major"),  fill = "transparent")

