
###################################################################
####Author:  Tatyana Doktorova
####Version: 1.0 12_12_2018
####Description: Disease selection
####            - filtering according to disease of interest
####            - this example is for HCC
####            - the lift can be adjusted appropriately in order to end up with similar number of output categories
#####################################################################

# Load
library("dplyr")
library(data.table)
library (VennDiagram)

##Set root location

Root_location <- "/home/dimiter/Desktop/Advance"


##Output directory
outDir <- paste(Root_location, "/OutputFiles/", sep='')

################################# Chemical X Gene data ##################
##Import the 3 overview files coming as a result of the ToxCast/CTD relationship code (https://github.com/DouglasConnect/advance).
####The files are available for download as well here: https://mega.nz/#F!zDBmiKyZ


TC <- read.table(paste(Root_location,"/OutputFiles/TC_ALLChem_Gene_protein_pathway_disease_associations.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', fill=TRUE, strip.white=TRUE, colClasses="character")

CTDTC <- read.table(paste(Root_location,"/OutputFiles/CTDTC_ALLChem_Gene_protein_pathway_disease_associations.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', fill=TRUE, strip.white=TRUE, colClasses="character")

CTD <- read.table(paste(Root_location,"/OutputFiles/CTDALL_ALLChem_Gene_protein_pathway_disease_associations.txt", sep=''), header=TRUE, sep='\t', comment.char='', quote='', fill=TRUE, strip.white=TRUE, colClasses="character")

##Filtering for HCC

TC1<-subset(TC, rhs=='D006528')
CTDTC1<-subset(CTDTC, rhs=='D006528')
CTD1<-subset(CTD, rhs=='D006528')


# Filtering for lift>50%
TC2<-transform (TC1, lift=as.numeric(lift))
TC3<- quantile(TC2$lift,  probs = c( 50)/100)
TC_perc<- data.table(TC2)
TC_perc<- TC_perc[lift > TC3]
write.table(TC_perc,file=file.path(outDir, "TC_perc.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

CTDTC2<-transform (CTDTC1, lift=as.numeric(lift))
CTDTC3<- quantile(CTDTC2$lift,  probs = c( 80)/100)
CTDTC_perc<- data.table(CTDTC2)
CTDTC_perc<- CTDTC_perc[lift > CTDTC3]
write.table(CTDTC_perc,file=file.path(outDir, "CTDTC_perc.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

CTD2<-transform (CTD1, lift=as.numeric(lift))
CTD3<- quantile(CTD2$lift,  probs = c( 95)/100)
CTD_perc<- data.table(CTD2)
CTD_perc<- CTD_perc[lift > CTD3]
write.table(CTD_perc,file=file.path(outDir, "CTD_perc.txt"), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)

#Unique Major/Minor only
TC_unique_major<-distinct(TC_perc, Major, .keep_all = TRUE)
TC_unique_minor<-distinct(TC_perc, Minor, .keep_all = TRUE)
CTDTC_unique_major<-distinct(CTDTC_perc, Major, .keep_all = TRUE)
CTDTC_unique_minor<-distinct(CTDTC_perc, Minor, .keep_all = TRUE)
CTD_unique_major<-distinct(CTDTC_perc, Major, .keep_all = TRUE)
CTD_unique_minor<-distinct(CTD_perc, Minor,  .keep_all = TRUE)


TC_unique_major_list<-distinct(TC_perc, Major )
TC_unique_minor_list<-distinct(TC_perc, Minor)
CTDTC_unique_major_list<-distinct(CTDTC_perc, Major)
CTDTC_unique_minor_list<-distinct(CTDTC_perc, Minor)
CTD_unique_major_list<-distinct(CTDTC_perc, Major)
CTD_unique_minor_list<-distinct(CTD_perc, Minor)

#Unique genes only

TC_unique<-distinct(TC_perc, lhs, .keep_all = TRUE)
CTDTC_unique<-distinct(CTDTC_perc, lhs, .keep_all = TRUE)
CTD_unique<-distinct(CTD_perc, lhs,  .keep_all = TRUE)

#Common genes/major/minor

Overlap_Genes<-intersect(TC_unique, CTDTC_unique, CTD_unique) 
Overlap_Major<-intersect(TC_unique_major_list, CTDTC_unique_major_list, CTD_unique_major_list) 
Overlap_Minor<-intersect(TC_unique_minor_list, CTDTC_unique_minor_list, CTD_unique_minor_list) 

