####################################################################
####Author: Noffisat Oki
####Version: 1.0 11-20-2017
####Description: Pubchem CID matching
####             Adds in Pubchem CIDS to the CTD and ToxCast data. 
####            This provides a common ID that can be used to match chemicals across both datasets 
####            as the pubchem SID is not sufficient for this purpose
####Notes: The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 15
####       should be used for this purpose.
####Potential issues:
#####################################################################

### Root folder location 
Root_location <- "C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

###reading mapping of all Pubchem IDs to CTD IDs#### the mapping file is obtained directly from PUBCHEM
ID_map_file_ctd <- read.delim(paste(Root_location, "/IntermediateData/pubchem2ctd_mapping.tsv", sep=''))
ids_ctd <- as.data.frame(ID_map_file_ctd[,c(3,4,5)]) 

###reading in mapping file of all Pubchem IDs to Toxcast GSIDs ###the mapping file is obtained directly from PUBCHEM
ID_map_file_toxcast <- read.delim(paste(Root_location, "/IntermediateData/pubchem2Toxcast_mapping.tsv", sep=''))
ids_toxcast <- as.data.frame(ID_map_file_toxcast[,c(3,4,5)]) 

#####Adding in all pubchem IDs to CTD Chem x Disease file
CTD_chemdisease_data <- read.csv(paste(Root_location,"/IntermediateData/CTD_chem_disease_wide_Dec172017download.csv",sep=''), header=TRUE, sep=",", fill=TRUE)
chem_disease_file_wt_all_IDs<-merge(ids_ctd, CTD_chemdisease_data, by.x="ctdid", by.y="CTD_ID")
write.table(chem_disease_file_wt_all_IDs,file=file.path(outDir, "CTD_chem_x_disease_wt_all_IDs.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)

#####Adding in all pubchem IDs to CTD Chem x Gene file (human gene set)
CTD_chemgene_data_human <- read.csv(paste(Root_location, "/IntermediateData/CTD_chem_gene_human_wide_Jan182018_download.csv", sep=''), header=TRUE, sep=",", fill=TRUE)
chem_gene_human_file_wt_all_IDs<-merge(ids_ctd, CTD_chemgene_data_human, by.x="ctdid", by.y="CTD_ID")
write.table(chem_gene_human_file_wt_all_IDs,file=file.path(outDir, "CTD_chem_x_gene_human_wt_all_IDs.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)

#####Adding in pubchem IDs to Toxcast Chem x Gene file for (human gene set)
Toxcast_chemgene_data_human <- read.csv(paste(Root_location,"/IntermediateData/ToxCastAssay_Hitcalls_chem_x_gene_wtPUBCHEM_human.csv",sep=''), header=TRUE, fill=TRUE)
Toxcast_chem_gene_file_wt_all_IDs_human<-merge(ids_toxcast, Toxcast_chemgene_data_human, by.x="sid", by.y="PUBCHEM_SID")
write.table(Toxcast_chem_gene_file_wt_all_IDs_human,file=file.path(outDir, "ToxCast_chem_x_gene_Hitcalls_wt_all_IDs_human.csv"), sep=',',col.names=TRUE, row.names=FALSE ,quote=FALSE)

