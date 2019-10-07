####PLEASE READ ALL NOTES and DESCRIPTIONS IN THIS HEADER BEFORE PROCEEDING TO THE CODE
##########################
####Author: Noffisat Oki
####Version: 1.0 02-02-2018
####Description: This script merges the mappings from genes to proteins and pathways between Uniprot and Reactome 
####             using the annotation information of the genes common to both the Toxcast and CTD datasets.
####
####Notes: This file is essentially a commented and cleaned up version of 'Gene_protein_Reactome_mapping.R'
####       An earlier version of this code that did not include the pathways (code after line 39) 
####       is 'Gene_protein_mapping.R'. 
####Potential issues:

#############################
##The root location folder should be specified here. 
#This is the main folder where any subfolders from which data will be read from or written to are located

Root_location <- "C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project/ADVANCE"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')
outDir2 <- paste(Root_location, "/OutputFiles/", sep='')

##Read in the Toxcast annotation data
Toxcast_annotation <- read.table(paste(Root_location,"/InputFiles/Toxcast_data_Oct2015_release/Assay_Information_Oct_2015/Assay_Summary_151020.csv", sep=''), header=TRUE, sep=',', comment.char='', quote="\"'", na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##Selelct only for the human gene data
Toxcast_annotation_human <- Toxcast_annotation[which(Toxcast_annotation$organism== "human"),]
##excluding extraneous columns
human_IDs <- na.omit(Toxcast_annotation_human[,c(55,59)])
##ensuring no duplicate genes
TC_human_IDs <- unique(human_IDs[which(human_IDs$intended_target_entrez_gene_id!="NA"),])
TC_human_IDs[TC_human_IDs=="HLA-DRA"] <- "HLA.DRA"


##Read in CTD gene annotation
################################# Chemical X Gene data ##################
##CTD chemical by gene data
CTD_chemgene_data.long <- read.table(paste(Root_location,"/InputFiles/CTD_chem_gene_ixns_Jan182018.tsv", sep=''), header=FALSE, sep='\t', skip=28, comment.char='', quote='', fill=TRUE, strip.white=TRUE, colClasses="character")
colnames(CTD_chemgene_data.long)<-  c("ChemicalName", "CTD_ID", "CASRN", "GeneSymbol", "GeneID", "GeneForm", "Organism", "OrganismID", "Interaction", "InteractionActions", "PubmedIDs")
##selecting human genes only
CTD_chemgene_data_human.long <- CTD_chemgene_data.long[which(CTD_chemgene_data.long$Organism=="Homo sapiens"),]
CTD_human_IDs <- unique(CTD_chemgene_data_human.long[,c(4,5)])



###Reading in the rules from the FIM analysis
######ToxCast
TC_FIM_analysis_genes <- read.table(paste(Root_location,"/OutputFiles/workflow4_Toxcast_new_stdlifts_wtDiseaseMap.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging the gene annotation data to the rules 
TC_FIM_genes_for_uniprot_mapping <- (merge(TC_human_IDs, TC_FIM_analysis_genes, by.x="intended_target_gene_symbol", by.y="lhs" ))
## keeping only required columns (gene symbol and entrez gene id)
TC_FIM_genes_for_uniprot_mapping <- unique(TC_FIM_genes_for_uniprot_mapping[,c(1,2)])
write.table(TC_FIM_genes_for_uniprot_mapping,file=file.path(outDir, "TC_FIMresults_genes_and_IDs.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
##The output directly above is used to get the mapping results 'Uniprot_gene_map_1_18_2018.tab' in the ../InputData/Uniprot directory


###Reading in the rules from the FIM analysis
######CTD (TC GENES)
CTDTC_FIM_analysis_genes <- read.table(paste(Root_location,"/OutputFiles/workflow4_TCgenes_in_CTD_new_stdlifts_wtDiseaseMap.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging the gene annotation data to the rules 
CTDTC_FIM_genes_for_uniprot_mapping <- (merge(CTD_human_IDs, CTDTC_FIM_analysis_genes, by.x="GeneSymbol", by.y="lhs" ))
## keeping only required columns (gene symbol and entrez gene id)
CTDTC_FIM_genes_for_uniprot_mapping <- unique(CTDTC_FIM_genes_for_uniprot_mapping[,c(1,2)])
write.table(CTDTC_FIM_genes_for_uniprot_mapping,file=file.path(outDir, "CTD_TC_FIMresults_genes_and_IDs.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
##The output directly above is used to get the mapping results 'Uniprot_gene_map_1_18_2018.tab' in the ../InputData/Uniprot directory

###Reading in the rules from the FIM analysis
######CTD ALL
CTDALL_FIM_analysis_genes <- read.table(paste(Root_location,"/OutputFiles/workflow4_allCTD_new_stdlifts_wtDiseaseMap.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging the gene annotation data to the rules 
CTDALL_FIM_genes_for_uniprot_mapping <- (merge(CTD_human_IDs, CTDALL_FIM_analysis_genes, by.x="GeneSymbol", by.y="lhs" ))
## keeping only required columns (gene symbol and entrez gene id)
CTDALL_FIM_genes_for_uniprot_mapping <- unique(CTDALL_FIM_genes_for_uniprot_mapping[,c(1,2)])
write.table(CTDALL_FIM_genes_for_uniprot_mapping,file=file.path(outDir, "CTD_ALL_FIMresults_genes_and_IDs.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
##The output directly above is used to get the mapping results 'Uniprot_gene_map_1_18_2018.tab' in the ../InputData/Uniprot directory



##Reading in the Uniprot (entrez gene id to uniprot protein ID) and Reactome (entrez gene id to protein identifier) mappings
Uniprot_gene_map_data <- read.table(paste(Root_location,"/IntermediateData/uniprot_gene_FIM_map_1_18_2018.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
Reactome_mapping_results_data <- read.table(paste(Root_location,"/IntermediateData/Reactome_analysis_mapping_FIM_TC_CTD_1_19_2018b.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging the uniprot and reactome maps
merged_gene_protein_human<- unique(merge(Uniprot_gene_map_data, Reactome_mapping_results_data, by.x="yourlist.M20180118A7434721E10EE6586998A056CCD0537E28D2D0H", by.y="Submitted.identifier"))
##merging the uniprot and reactome IDs to the FIM rules and annotation data
###Toxcast
TC_merged_with_GeneID <- merge(TC_FIM_genes_for_uniprot_mapping,merged_gene_protein_human, by.x= "intended_target_entrez_gene_id", by.y="yourlist.M20180118A7434721E10EE6586998A056CCD0537E28D2D0H") 
colnames(TC_merged_with_GeneID) <- c("intended_target_entrez_gene_id","intended_target_gene_symbol","Entry.Uniprot","Found.identifier.Reactome")#, "Resource")
write.table(TC_merged_with_GeneID,file=file.path(outDir, "Gene_protein_maps_TC_FIMresults.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
###CTD (TC GENES)
CTDTC_merged_with_GeneID <- merge(CTDTC_FIM_genes_for_uniprot_mapping,merged_gene_protein_human, by.x= "GeneID", by.y="yourlist.M20180118A7434721E10EE6586998A056CCD0537E28D2D0H") 
colnames(CTDTC_merged_with_GeneID) <- c("intended_target_entrez_gene_id","intended_target_gene_symbol","Entry.Uniprot","Found.identifier.Reactome")#, "Resource")
write.table(CTDTC_merged_with_GeneID,file=file.path(outDir, "Gene_protein_maps_CTDTC_FIMresults.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
###CTD ALL
CTDALL_merged_with_GeneID <- merge(CTDALL_FIM_genes_for_uniprot_mapping,merged_gene_protein_human, by.x= "GeneID", by.y="yourlist.M20180118A7434721E10EE6586998A056CCD0537E28D2D0H") 
colnames(CTDALL_merged_with_GeneID) <- c("intended_target_entrez_gene_id","intended_target_gene_symbol","Entry.Uniprot","Found.identifier.Reactome")#, "Resource")
write.table(CTDALL_merged_with_GeneID,file=file.path(outDir, "Gene_protein_maps_CTDALL_FIMresults.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)


##Reading in the full reactome to uniprot pathway data. This is output from the 'Merge_analysis_1.R' (ReactomeClassv2.R) code. See the script
##for more information.
ReactomeToUniprot <- read.table(paste(Root_location,"/IntermediateData/ReactomePathways2UniProt.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote="", na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##Selecting only for human pathways
ReactomeToUniprot_human <- ReactomeToUniprot[which(ReactomeToUniprot$Species=="Homo sapiens"),]
##Re-reading in the file of the FIM rules and annotation with protein identifier info created above
##TOXCAST
TC_Gene_protein <- read.table(paste(Root_location,"/IntermediateData/Gene_protein_maps_TC_FIMresults.csv", sep=''),header=TRUE, sep=',', comment.char='', quote="", na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging this with the merged FIM annonation and protein file from above
TC_merged_gene_protein_reactome<-(merge(TC_Gene_protein, ReactomeToUniprot_human, by.x="Found.identifier.Reactome", by.y="UniProtID"))
write.table(TC_merged_gene_protein_reactome,file=file.path(outDir, "TC_Gene_protein_pathway_maps.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
##CTD (TC GENES)
CTDTC_Gene_protein <- read.table(paste(Root_location,"/IntermediateData/Gene_protein_maps_CTDTC_FIMresults.csv", sep=''),header=TRUE, sep=',', comment.char='', quote="", na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging this with the merged FIM annonation and protein file from above
CTDTC_merged_gene_protein_reactome<-(merge(CTDTC_Gene_protein, ReactomeToUniprot_human, by.x="Found.identifier.Reactome", by.y="UniProtID"))
write.table(CTDTC_merged_gene_protein_reactome,file=file.path(outDir, "CTDTC_Gene_protein_pathway_maps.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
##CTD ALL
CTDALL_Gene_protein <- read.table(paste(Root_location,"/IntermediateData/Gene_protein_maps_CTDALL_FIMresults.csv", sep=''),header=TRUE, sep=',', comment.char='', quote="", na.strings = "",fill=TRUE, strip.white=TRUE, colClasses="character")
##merging this with the merged FIM annonation and protein file from above
CTDALL_merged_gene_protein_reactome<-(merge(CTDALL_Gene_protein, ReactomeToUniprot_human, by.x="Found.identifier.Reactome", by.y="UniProtID"))
write.table(CTDALL_merged_gene_protein_reactome,file=file.path(outDir, "CTDALL_Gene_protein_pathway_maps.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)



########Merging with disease information#########
TC_gene_protein_pathway_disease <- unique(merge(TC_FIM_analysis_genes, TC_merged_gene_protein_reactome, by.x="lhs", by.y="intended_target_gene_symbol"))
write.table(TC_gene_protein_pathway_disease, file=file.path(outDir2, "TC_Gene_protein_pathway_disease_associations.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

CTDTC_gene_protein_pathway_disease <- unique(merge(CTDTC_FIM_analysis_genes, CTDTC_merged_gene_protein_reactome, by.x="lhs", by.y="intended_target_gene_symbol"))
write.table(CTDTC_gene_protein_pathway_disease, file=file.path(outDir2, "CTDTC_Gene_protein_pathway_disease_associations.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

CTDALL_gene_protein_pathway_disease <- unique(merge(CTDALL_FIM_analysis_genes, CTDALL_merged_gene_protein_reactome, by.x="lhs", by.y="intended_target_gene_symbol"))
write.table(CTDALL_gene_protein_pathway_disease, file=file.path(outDir2, "CTDALL_Gene_protein_pathway_disease_associations.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
