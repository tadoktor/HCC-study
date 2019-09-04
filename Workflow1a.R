####Author: Noffisat Oki
####Version: 1.0 11-20-2017
####Description: Reformatting CTD data from long to wide format
####             Converts CTD chemical-gene and chemical-disease data from long to wide 
####             format accepted by the 'arules' package. 
####Notes: Current version selects only human (homo sapiens) data from the CTD dataset. 
####       Script may be modified for analysis of other desired species present in the CTD data.
####       The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 19
####       should be used for this purpose.
####Potential issues:no known issues
###
###Removed extraneous code after line 43
#####################################################################
library(plyr)
library(data.table)
library(reshape2)

##The root location folder should be specified here. 
#This is the main folder where any subfolders from which data will be read from or written to are located
Root_location <- "C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

################################# Chemical X Gene data ##################
##CTD chemical by gene data
CTD_chemgene_data.long <- read.table(paste(Root_location,"/InputFiles/CTD_chem_gene_ixns_Jan182018.tsv", sep=''), header=FALSE, sep='\t', skip=28, comment.char='', quote='', fill=TRUE, strip.white=TRUE, colClasses="character")
colnames(CTD_chemgene_data.long)<-  c("ChemicalName", "CTD_ID", "CASRN", "GeneSymbol", "GeneID", "GeneForm", "Organism", "OrganismID", "Interaction", "InteractionActions", "PubmedIDs")
##selecting human genes only
CTD_chemgene_data_human.long <- CTD_chemgene_data.long[which(CTD_chemgene_data.long$Organism=="Homo sapiens"),]
##reformatting the data
CTD_chemgene_data_human.wide <- dcast(CTD_chemgene_data_human.long, CTD_ID ~ GeneSymbol , value.var="GeneID")
write.table(CTD_chemgene_data_human.wide,file=file.path(outDir, "CTD_chem_gene_human_wide_Jan182018_download.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)

######################  Chemical X Disease  ###############
##CTD chemical by disease data
#CTD_chemdisease_data.long <- read.table("/wikihomedir/noki/cpAOP_paper/InputData/CTD_chemicals_diseases_2_27_2015.txt", header=FALSE, sep='\t', skip=28, comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
CTD_chemdisease_data.long <- read.table(paste(Root_location,"/InputFiles/CTD_chemicals_diseases_Dec192017.tsv", sep=''), header=FALSE, sep='\t', skip=28, comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
colnames(CTD_chemdisease_data.long)<-  c("ChemicalName", "CTD_ID", "CASRN", "DiseaseName", "DiseaseID", "DirectEvidence", "InferenceGeneSymbol", "InferenceScore", "OmimIDs", "PubmedIDs")
##reformatting the data
CTD_chemdisease_data.wide <- dcast(CTD_chemdisease_data.long, CTD_ID ~ DiseaseID ,value.var="DiseaseName")
write.table(CTD_chemdisease_data.wide,file=file.path(outDir, "CTD_chem_disease_wide_Dec172017download.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)

