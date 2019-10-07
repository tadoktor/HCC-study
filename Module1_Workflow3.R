####################################################################
####Author: Noffisat Oki
####Version: 1.0 11-20-2017
####Description: Gene matching
####             Matches genes between ToxCast and CTD data.
####Notes: Current version analyzes only human data from both the CTD and ToxCast datasets. 
####       Script may be modified for analysis of other desired species present in the ToxCast data.
####       The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 15
####       should be used for this purpose.
####Potential issues:
#####################################################################

### Root folder location 
Root_location <- "C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

####selecting and Matching genes from CTD to ToxCast ##########################
###ToxCast Human gene targets
ToxCast_chemgene_wt_PubchemID_human <- read.csv(paste(Root_location,"/IntermediateData/ToxCast_chem_x_gene_Hitcalls_wt_all_IDs_human.csv", sep=''), header=TRUE, sep=",", fill=TRUE)
genes_in_toxcast_human <- colnames(ToxCast_chemgene_wt_PubchemID_human)
genes_in_toxcast_human <- genes_in_toxcast_human[c(-1,-2,-3)]
##CTD Human gene set
CTD_chemgene_wt_PubchemID_human <- read.csv(paste(Root_location, "/IntermediateData/CTD_chem_x_gene_human_wt_all_IDs.csv",sep=''), header=TRUE, sep=",", fill=TRUE)
##matching genes between datasets
myvars2 <- which(colnames(unique(CTD_chemgene_wt_PubchemID_human)) %in% genes_in_toxcast_human)
match_CTD_chemgene_data_human <- (CTD_chemgene_wt_PubchemID_human[,c(1:3,myvars2)])
write.table(match_CTD_chemgene_data_human,file=file.path(outDir, "CTD_chem_x_gene_human_wt_all_IDs_humanTCgenes.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
