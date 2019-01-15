####################################################################
####Author: Noffisat Oki, Tatyana Doktorova
####Version: 1.0 08-28-2018
####Description: GTX Subset script:
####             Subsets chemical-gene data to filter for only genes active
####             in at least 3 chemicals 
###               Reannotation to Ensemb Identifier and visualizations
####             
####Notes: The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 20
####       should be used for this purpose.
####Potential issues:
#####################################################################

library(plyr)
library(data.table)
library(reshape2)

##The root location folder should be specified here. 
#This is the main folder where any subfolders from which data will be read from or written to are located
Root_location <- "C:\\Users\\tatyana\\Documents\\Projects\\Advance project\\Advance"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

################################# Chemical X Gene data ##################
GTX_chemgene_data.long <- read.csv(paste(Root_location,"/InputFiles/GTX_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
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
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

mapping_GTX <- getBM(attributes = c("affy_hg_u133_plus_2",
                                    "ensembl_gene_id", "hgnc_symbol"), filters = "affy_hg_u133_plus_2",
                     values = Probes, mart = ensembl)


GTX_GL_final <- mapping_GTX[!duplicated(mapping_GTX[,3]),] 


################################# NGTX Chemical X Gene data ##################
NGTX_chemgene_data.long <- read.csv(paste(Root_location,"/InputFiles/NGTX_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
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




mapping_NGTX <- getBM(attributes = c("affy_hg_u133_plus_2",
                                    "ensembl_gene_id", "hgnc_symbol"), filters = "affy_hg_u133_plus_2",
                     values = NGTX_Probes, mart = ensembl)


NGTX_GL_final <- mapping_NGTX[!duplicated(mapping_NGTX[,3]),] 



################################# NC Chemical X Gene data ##################
NC_chemgene_data.long <- read.csv(paste(Root_location,"/InputFiles/NC_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
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




mapping_NC <- getBM(attributes = c("affy_hg_u133_plus_2",
                                     "ensembl_gene_id", "hgnc_symbol"), filters = "affy_hg_u133_plus_2",
                      values = NC_Probes, mart = ensembl)


NC_GL_final <- mapping_NC[!duplicated(mapping_NC[,3]),] 





GTX_GL_genes<-unique(subset(GTX_GL_final, select=c("hgnc_symbol")))
NGTX_GL_genes<-unique(subset(NGTX_GL_final, select=c("hgnc_symbol")))
NC_GL_genes<-unique(subset(NC_GL_final, select=c("hgnc_symbol")))


Genes_common<-intersect(GTX_GL_genes, NGTX_GL_genes,NC_GL_genes)


NGTX_specific_1<-setdiff(NGTX_GL_genes, GTX_GL_genes)
NGTX_specific_2<-setdiff(NGTX_GL_genes, NC_GL_genes)
NGTX_specific<-merge(NGTX_specific_1, NGTX_specific_2)

NGTX_specific<-data.frame(NGTX_specific)


HCC_Toxcast_CTD<-read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/HCC_Toxcast_CTD.csv", header=TRUE)
Unique_genes = HCC_Toxcast_CTD %>% distinct(lhs)
colnames(Unique_genes) <- c("hgnc_symbol")
Unique_genes<-data.frame(Unique_genes)

Genes_TGGATEs_CTD_Toxcast<-intersect(NGTX_specific, Unique_genes)


NGTX_induced_HCC_overlap<-merge(Genes_TGGATEs_CTD_Toxcast, HCC_Toxcast_CTD, by.x="hgnc_symbol", by.y="lhs")

NGTX<-rbind(NGTX_induced_HCC,NGTX_induced_HCC_overlap)

collapsibleTree(df=NGTX, c("hgnc_symbol", "Major","Minor","Parent.Name"),  fill = "transparent")

collapsibleTree(df=NGTX_induced_HCC_overlap, c("Major", "Minor","Parent.Name", "hgnc_symbol"),  fill = "lightgreen")


collapsibleTree(df=NGTX_induced_HCC_overlap, c("hgnc_symbol", "Parent.Name"),  fill = "transparent")

collapsibleTree(df=NGTX_induced_HCC_overlap, c("hgnc_symbol", "Major"),  fill = "transparent")
