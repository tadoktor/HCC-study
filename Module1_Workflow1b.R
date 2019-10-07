####################################################################
####Author: Noffisat Oki
####Version: 1.0 11-20-2017
####Description: ToxCast Mapping script:
####             Maps ToxCast chemical-assay hitcall data to chemical-gene 
####             
####Notes: Current version selects only human assay data from the ToxCast dataset. 
####       Script may be modified for analysis of other desired species present in the ToxCast data.
####       The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 16
####       should be used for this purpose.
####Potential issues:
#####################################################################

library(gdata)
### Root folder location 
Root_location <- "C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

##Read in Toxcast assay hitcall matrix file
Assay_data <- read.csv(paste(Root_location,"/InputFiles/Toxcast_data_Oct2015_release/INVITRODB_V2_SUMMARY/hitc_Matrix_151020.csv", sep=''), header=TRUE, sep=",", row.names=1,fill=TRUE)
##0/1 discretization (converting NAs and -1s to 0)
Assay_data_noNA <- Assay_data
Assay_data_noNA[is.na(Assay_data_noNA)] <- 0
Assay_data_noNA[Assay_data_noNA==-1] <- 0

####Subsetting for human assay data using the annotation study design info file
Annotation_file <- read.csv(paste(Root_location,"/InputFiles/Toxcast_data_Oct2015_release/Assay_Information_Oct_2015/Assay_Summary_151020.csv", sep=''), header=TRUE, sep=",", fill=TRUE)
human_assay_names <- Annotation_file[which(Annotation_file$organism=='human'),]
human_assay_names <- as.character(human_assay_names$assay_component_endpoint_name)
human_assay_data <- Assay_data_noNA[,which(c(colnames(Assay_data_noNA))%in% (human_assay_names))]

##Mapping the assay names to their gene targets
AssayGene_MapValues <- cbind(as.character(Annotation_file$assay_component_endpoint_name), toupper(as.character(Annotation_file$intended_target_gene_symbol)))
newAssay_colnames_human <- mapvalues(colnames(human_assay_data), from=c(AssayGene_MapValues[,1]), to=c(AssayGene_MapValues[,2]))
colnames(human_assay_data)<- newAssay_colnames_human

##aggregating the data Reformatting the table so that each gene is represented by one column
reshaping_assaydata_human <- t(human_assay_data)
reshaping_assaydata_agg_human <- sapply(by(reshaping_assaydata_human, rownames(reshaping_assaydata_human),colSums), identity)
##converting back to (0/1)
reshaping_assaydata_agg_human[reshaping_assaydata_agg_human>1] <- 1

ToxCast_chem_gene_human <- as.data.frame(reshaping_assaydata_agg_human)
write.csv(ToxCast_chem_gene_human,file=file.path(outDir, "ToxCastAssay_Hitcalls_chem_x_gene_human.csv"), col.names=TRUE, row.names=TRUE,quote=FALSE)

###Reading chemical ID mapping information
#assay_chemical_summary_file <- read.csv(paste(Root_location, "/InputFiles/Chemical_Summary_151020.csv", sep=''), header=TRUE, sep=",", fill=TRUE)
assay_chemical_summary_file <- read.csv(paste(Root_location, "/InputFiles/Toxcast_data_Oct2015_release/INVITRODB_V2_SUMMARY/Chemical_Summary_151020.csv", sep=''), header=TRUE, sep=",", fill=TRUE)
chemical_summary_file2 <- read_excel(paste(Root_location, "/InputFiles/Toxcast_data_Oct2015_release/Chemical_information/TOX21IDs_v4b_23Oct2014_QCdetails.xlsx", sep=''), sheet=1)


##ADD in DSSTox GSIDs
Casrn_codes_human <- rownames(ToxCast_chem_gene_human)
DSSTox_IDs_human <- as.data.frame(mapvalues(Casrn_codes_human, from=as.character((assay_chemical_summary_file$code)), to=as.character((assay_chemical_summary_file$chid))))
colnames(DSSTox_IDs_human)<- "DSSTox_GSID"
ToxCast_chem_gene_human2 <- cbind(DSSTox_IDs_human,ToxCast_chem_gene_human)
write.table(ToxCast_chem_gene_human2,file=file.path(outDir, "ToxCastAssay_Hitcalls_chem_x_gene_wtDSSTOX_human.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

#####Adding in pubchem_SIDs
chemIDs_human <- as.character(ToxCast_chem_gene_human2$DSSTox_GSID)
Pubchem_SIDs_human <- as.data.frame(mapvalues(chemIDs_human, from=c(chemical_summary_file2$ToxCast_chid), to=c(chemical_summary_file2$PUBCHEM_SID)))
colnames(Pubchem_SIDs_human) <- "PUBCHEM_SID"
chem_gene_file_wt_pubchemSID_human <- cbind(Pubchem_SIDs_human, ToxCast_chem_gene_human2)
write.csv(chem_gene_file_wt_pubchemSID_human,file=file.path(outDir, "ToxCastAssay_Hitcalls_chem_x_gene_wtPUBCHEM_human.csv"), col.names=TRUE, row.names=TRUE,quote=FALSE)

