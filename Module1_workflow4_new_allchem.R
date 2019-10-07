####################################################################
####Author: Noffisat Oki
####Edits: 
####Version: 2.0 11-22-2017
####Description: Association rule analysis workflow.
####Notes: Current version analyzes only human data from both the CTD and ToxCast datasets. 
####       Modifications have been made to the arules_function so that it can run on data of "any" size, to avoid
####       having to run the function on multiple subsets of a large dataframe because of memory limitations.  Other modifications
####       include the addition of two different standarized lift values.
####       The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. 
####Potential issues:
#####################################################################


library(gdata)
library(arules)
#setwd("C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project")

### Root folder location 
Root_location <- "C:/Users/nofis/Documents/DouglasConnect/ADVANCE_Project/ADVANCE"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')


Gene_data_Toxcast_human1 <- read.csv(paste(Root_location,"/IntermediateData/ToxCast_chem_x_gene_Hitcalls_wt_all_IDs_human.csv", sep=''),sep=',', header=TRUE, comment.char='', strip.white=TRUE)
Idcols <- Gene_data_Toxcast_human1[,c(1:5)]
Gene_data_Toxcast_humanActive <-  Gene_data_Toxcast_human1
Gene_data_Toxcast_human_unique <- unique(Gene_data_Toxcast_human1[,c(-1,-2,-3,-4)])
Gene_data_Toxcast_human <- Gene_data_Toxcast_human1[,c(-1,-3,-4,-5)]
###Remove chemicals with no CID (CID=NA) as the gene and disease datasets are merged on CID
Gene_data_Toxcast_human <- na.omit(Gene_data_Toxcast_human)
ids <- Gene_data_Toxcast_human[,1]
##making sure all values are discretized either 1 or 0
Gene_data_Toxcast_human[Gene_data_Toxcast_human>=1] <- 1
Gene_data_Toxcast_human[Gene_data_Toxcast_human<1] <- 0
Gene_data_Toxcast_human <- Gene_data_Toxcast_human[,c(-1)]
Gene_data_Toxcast_human <- cbind(ids, Gene_data_Toxcast_human)
colnames(Gene_data_Toxcast_human)[1]<-"cid"


####To print out the active genes and chemicals matrix for Toxcast set
###Remove chemicals with no CID (CID=NA) as the gene and disease datasets are merged on CID
Gene_data_Toxcast_humanActive <- na.omit(Gene_data_Toxcast_humanActive)
ids2 <- Gene_data_Toxcast_humanActive[,1:5]
Gene_data_Toxcast_humanActive2 <- Gene_data_Toxcast_humanActive[,c(-1,-2,-3,-4,-5)]
##making sure all values are discretized either 1 or 0
Gene_data_Toxcast_humanActive2[Gene_data_Toxcast_humanActive2>=1] <- 1
Gene_data_Toxcast_humanActive2[Gene_data_Toxcast_humanActive2<1] <- 0
rownames(Gene_data_Toxcast_humanActive2) <- ids2[,1]
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_dataset<-Gene_data_Toxcast_humanActive2[,!apply(Gene_data_Toxcast_humanActive2,2,function(x) sum(abs(x), na.rm=TRUE) < 1)]
#remove rows which are all '0' indicating no effect for a chemical
Active_dataset<-Active_dataset[!apply(Active_dataset,1,function(x) sum(abs(x), na.rm=TRUE)< 1),]
sid <- as.matrix(rownames(Active_dataset))
Active_dataset <- cbind(sid,Active_dataset)
Active_dataset <- merge(ids2,Active_dataset, by.x='sid', by.y='sid')
dataset <- "/OutputFiles/Toxcast_ALLchems_Active_dataset.txt"
write.table(Active_dataset,paste(Root_location,dataset,sep=''),sep="\t",row.names=F)




##TC genes in CTD
Gene_data_CTD_TC_human1 <-read.csv(paste(Root_location,"/IntermediateData/CTD_chem_x_gene_human_wt_all_IDs_humanTCgenes.csv", sep=''),sep=",",header=T,comment.char='',strip.white=T)
Gene_data_CTD_TC_human <- na.omit(Gene_data_CTD_TC_human1)
Gene_data_CTD_TC_humanActive <- Gene_data_CTD_TC_human
Gene_data_CTD_TC_human_unique <- unique(Gene_data_CTD_TC_human[,c(-1,-2)])
Gene_data_CTD_TC_human <- unique(Gene_data_CTD_TC_human[,c(-1,-2)])
cid <- Gene_data_CTD_TC_human[,1]
Gene_data_CTD_TC_human <- Gene_data_CTD_TC_human[,c(-1)]
##making sure all values are discretized either 1 or 0
Gene_data_CTD_TC_human[Gene_data_CTD_TC_human>=1] <- 1
Gene_data_CTD_TC_human[Gene_data_CTD_TC_human<1] <- 0
Gene_data_CTD_TC_human <- cbind(cid, Gene_data_CTD_TC_human)



####Print out the active genes and chemicals matrix for CTD Toxcast set
###Remove chemicals with no CID (CID=NA) as the gene and disease datasets are merged on CID
Gene_data_CTD_TC_humanActive <- na.omit(Gene_data_CTD_TC_humanActive)
ids3 <- Gene_data_CTD_TC_humanActive[,1:3]
Gene_data_CTD_TC_humanActive2 <- Gene_data_CTD_TC_humanActive[,c(-1,-2,-3)]
##making sure all values are discretized either 1 or 0
Gene_data_CTD_TC_humanActive2[Gene_data_CTD_TC_humanActive2>=1] <- 1
Gene_data_CTD_TC_humanActive2[Gene_data_CTD_TC_humanActive2<1] <- 0
uniqueIDs <- as.data.frame(matrix(c(1:nrow(ids3)),nrow=nrow(ids3),ncol=1) )
uniqueIDs <- cbind (uniqueIDs,ids3)
rownames(Gene_data_CTD_TC_humanActive2) <- uniqueIDs[,1]
#remove columns which are all '0' indicating no effect for an assay, gene or disease
CTDTC_Active_dataset<- Gene_data_CTD_TC_humanActive2[,!apply(Gene_data_CTD_TC_humanActive2,2,function(x) sum(abs(x), na.rm=TRUE) < 1)]
#remove rows which are all '0' indicating no effect for a chemical
CTDTC_Active_dataset<-CTDTC_Active_dataset[!apply(CTDTC_Active_dataset,1,function(x) sum(abs(x), na.rm=TRUE)< 1),]
uid <- as.matrix(rownames(CTDTC_Active_dataset))
CTDTC_Active_dataset <- cbind(uid,CTDTC_Active_dataset)
colnames(uniqueIDs)<-c("UID","CTDID","PubChem_SID", "PubChem_CID")
CTDTC_Active_dataset <- merge(uniqueIDs, CTDTC_Active_dataset, by.x='UID', by.y='uid')
dataset <- "/OutputFiles/CTDTC_ALLchems_Active_dataset.txt"
write.table(CTDTC_Active_dataset,paste(Root_location,dataset,sep=''),sep="\t",row.names=F)



##All genes in CTD
CTD_all1 <-read.csv(paste(Root_location,"/IntermediateData/CTD_chem_x_gene_human_wt_all_IDs.csv",sep=''),sep=",",header=T,comment.char='', strip.white=T)
CTD_all<-na.omit(CTD_all1)
Gene_data_CTDALL_humanActive <- CTD_all
CTD_all_unique <- unique(CTD_all[,c(-1,-2)])
CTD_all<- CTD_all[,c(-1,-2)]
cid <- CTD_all[,1]
CTD_all <- CTD_all[,c(-1)]
##making sure all values are discretized either 1 or 0
CTD_all[CTD_all>=1] <- 1
CTD_all[CTD_all<1] <- 0
CTD_all <- cbind(cid, CTD_all)



####Print out the active genes and chemicals matrix for CTDALL set
###Remove chemicals with no CID (CID=NA) as the gene and disease datasets are merged on CID
Gene_data_CTDALL_humanActive <- na.omit(Gene_data_CTDALL_humanActive)
ids4 <- Gene_data_CTDALL_humanActive[,1:3]
Gene_data_CTDALL_humanActive2 <- Gene_data_CTDALL_humanActive[,c(-1,-2,-3)]
##making sure all values are discretized either 1 or 0
Gene_data_CTDALL_humanActive2[Gene_data_CTDALL_humanActive2>=1] <- 1
Gene_data_CTDALL_humanActive2[Gene_data_CTDALL_humanActive2<1] <- 0
uniqueIDs2 <- as.data.frame(matrix(c(1:nrow(ids4)),nrow=nrow(ids4),ncol=1) )
uniqueIDs2 <- cbind (uniqueIDs2,ids4)
rownames(Gene_data_CTDALL_humanActive2) <- uniqueIDs2[,1]
#rownames(Gene_data_CTDALL_humanActive2) <- ids4[,2]
#colnames(Gene_data_Toxcast_human)[1]<-"sid"
#remove columns which are all '0' indicating no effect for an assay, gene or disease
CTDALL_Active_dataset<-Gene_data_CTDALL_humanActive2[,!apply(Gene_data_CTDALL_humanActive2,2,function(x) sum(abs(x), na.rm=TRUE) < 1)]
#remove rows which are all '0' indicating no effect for a chemical
CTDALL_Active_dataset<-CTDALL_Active_dataset[!apply(CTDALL_Active_dataset,1,function(x) sum(abs(x), na.rm=TRUE)< 1),]
uid <- as.matrix(rownames(CTDALL_Active_dataset))
CTDALL_Active_dataset <- cbind(uid,CTDALL_Active_dataset)
colnames(uniqueIDs2)<-c("UID","CTDID","PubChem_SID", "PubChem_CID")
CTDALL_Active_dataset <- merge(uniqueIDs2,CTDALL_Active_dataset, by.x='UID', by.y='uid')
dataset <- "/OutputFiles/CTDALL_ALLchem_Active_dataset.txt"
write.table(CTDALL_Active_dataset,paste(Root_location,dataset,sep=''),sep="\t",row.names=F)



##All diseases in CTD
Disease_data_CTD1 <- read.csv(paste(Root_location,"/IntermediateData/CTD_chem_x_disease_wt_all_IDs.csv", sep=''), sep=',', header=TRUE, comment.char='', strip.white=TRUE)
Disease_data_CTD <- na.omit(Disease_data_CTD1)
Disease_data_CTD_Unique <- unique(Disease_data_CTD[,c(-1,-2)])
Disease_data_CTD <- Disease_data_CTD[,c(-1,-2)]
cids <- Disease_data_CTD[,1]
Disease_data_CTD <- Disease_data_CTD[,c(-1)]
##making sure all values are discretized either 1 or 0
Disease_data_CTD[Disease_data_CTD>=1] <- 1
Disease_data_CTD[Disease_data_CTD<1] <- 0
Disease_data_CTD <- cbind(cids, Disease_data_CTD)
colnames(Disease_data_CTD)[1]<-"cid"



#dataset2 <- "/OutputFiles/CTD_Disease_HCCcmpds_dataset.txt"
#write.table(Disease_data_CTD,paste(Root_location,dataset2,sep=''),sep="\t",row.names=F)


arules_function <- function(AssayData,PhenoData,sup,minactivity,logfile){ #}, datamerge){
  
  #AssayData<- Gene_data_Toxcast_human
  #PhenoData <-Disease_data_CTD
  #sup <-10
  #logfile <- logfile3
  
  assays<-colnames(AssayData)
  numgenes <- length(assays)-1
  pheno<-colnames(PhenoData) 
  numpheno <- length(pheno)-1
  dataDesc<-(merge(AssayData, PhenoData, by.x="cid", by.y="cid"))
  dataDesc <- subset(dataDesc, !duplicated(dataDesc[,1])) 
  
  
  
  #Must move the labels to the outside/as row names
  dataDesc2<-subset(dataDesc, select =-cid); row.names(dataDesc2)<-dataDesc$cid
  
  #to remove missing/uninformative datapoints.
  #remove COLUMNS that are missing all values
  dataBin<-dataDesc2[,!apply(is.na(dataDesc2),2,all)]
  #now remove rows that are missing all values
  #(can move to ANY if you are interested in a specific chemical)
  dataBin<-dataBin[!apply(is.na(dataBin),1,all),]
  
  dataBin[is.na(dataBin)]<-0
  
  #remove columns which are all '0' indicating no effect for an assay, gene or disease
  dataBinS<-dataBin[,!apply(dataBin,2,function(x) sum(abs(x), na.rm=TRUE) < minactivity)]
  #remove rows which are all '0' indicating no effect for a chemical
  dataBinS<-dataBinS[!apply(dataBinS,1,function(x) sum(abs(x), na.rm=TRUE)< minactivity),]
  
  #write.table(dataBinS,paste(Root_location,dataset,sep=''),sep="\t",row.names=F)
  
  
  ###Adding these four extra lines to eliminate the possibility of inf when calculating odds ratio
  combined_data <- colnames(dataBinS)
  genes<- which(assays%in%colnames(dataBinS))
  diseases<- which(pheno%in%colnames(dataBinS))
  tmp1 <- matrix(nrow=1, ncol=ncol(dataBinS), data=0)
  tmp1[,(length(genes)+1):ncol(dataBinS)] <- 1
  tmp2 <- matrix(nrow=1, ncol=ncol(dataBinS), data=0)
  tmp2[,1:length(genes)] <- 1
  tmp3 <- matrix(nrow=1, ncol=ncol(dataBinS), data=0)
  tmp4 <- matrix(nrow=1, ncol=ncol(dataBinS), data=1)
  colnames(tmp1)<- colnames(dataBinS)
  colnames(tmp2)<- colnames(dataBinS)
  colnames(tmp3)<- colnames(dataBinS)
  colnames(tmp4)<- colnames(dataBinS)
  dataBinS<-rbind(dataBinS, tmp1,tmp2,tmp3,tmp4)
  
  #Root_locationx <- "C:/Users/Noffisat/Documents/DouglasConnect/ADVANCE_Project/ADVANCE/OutputFiles/"
  #write.table(dataDesc,paste(Root_locationx,datamerge,sep=''),sep="\t",row.names=F)
  
  ###Writing QC info to logfile 
  sink(file=logfile, append=TRUE)
  date <- Sys.Date()
  time <- Sys.time()
  today <- format(date, format="%B %d %Y")
  print(today)
  print(time)
  print("dimensions of merged data")
  print(dim(dataDesc2))
  print("dimensions of QCed data")
  print(dim(dataBinS))
  print("Number of transactions in QCed data")
  print(nrow(dataBinS))
  combined_data <- colnames(dataBinS)
  genes<- which(assays%in%colnames(dataBinS))
  diseases<- which(pheno%in%colnames(dataBinS))
  print("number of QCed genes")
  print(length(genes))
  print("number of QCed diseases")
  print(length(diseases))  
  print("minimum support (minimum number of transactions)")
  print(sup)
  print("minimum support %")
  sup <- sup/nrow(dataBinS)
  
  
  #to simplify the rules, removing the level '0' and replacing it with NA
  #the effect will be to have just presence for a 'response' measure
  dataBinEffect<-dataBinS;
  print("dim of dataBinS")
  print(dim(dataBinS))
  
  #dataBinS<-as.data.frame(apply(dataBinS,2,as.factor))   ###is this line needed?? comment out for now
  print("checkpoint 0")
  
  dataBinEffect[dataBinEffect==0]<-NA;
  
  
  
  data_stdlifts<-dataBinEffect
  
  
  
  dataBinEffect<-as.data.frame(apply(dataBinEffect,2,as.factor))
  
  
  print("checkpoint 0a")
  #create the transaction object
  bin.data<-as(dataBinEffect, "transactions")
  print("checkpoint 0b")
  #bin.data2<-as(dataBinEffect2, "transactions")
  #bin.data3<-as(dataBinEffect3, "transactions")
  #bin.data4<-as(dataBinEffect4, "transactions")
  #bin.data5<-as(dataBinEffect5, "transactions")
  #bin.data6<-as(dataBinEffect6, "transactions")
  
  #this gets out the "labels" for the rules coming from each dataset
  assayItems<-intersect(colnames(bin.data), paste(assays,1, sep='='))
  phenoItems<-intersect(colnames(bin.data), paste(pheno,1, sep='='))
  
  #############################################
  #Support = #instance of your 'set'/ total # transactions
  #Confidence = support(x,y)/support(LHS)
  #Lift = support(x,y)/support(x)*support(y) (ratio of observed support to that if X and Y were independent)
  #setting the min/max len=2 will limit to 2 itemset rules
  print("Checkpoint 1a")
  rules<-apriori(bin.data, parameter =list(support = sup, confidence = 0.001, minlen=2, maxlen=2))
  print("Checkpoint 1")
  print("Number of rules before subset")
  print(dim(rules))
  rules.ap<-subset(rules, subset = lift >1)
  print("Checkpoint 2")
  print("Number of rules after subset")
  
  
  print(class(rules.ap))
  print(length(rules.ap))  ####there are close to 5.5 million rules in rules.ap when using the full CTD gene database
  
  rules_split<-list()
  
  ###################################################################
  if (length(rules.ap)<6000000){
    rules_split<-split(rules.ap,cut(seq(length(rules.ap)),5))  ### this cut function ensures that the remainder of rules.ap after the split is spread across the splits as evenly as possible 
  }
  
  else rules_split<-split(rules.ap,cut(seq(length(rules.ap)),10))
  ###################################################################   ^this section needs to be developed more, how to decide how many splits should be made?  
  #The number of splits doesn't affect the output so maybe just leave it as a large enough constant
  #so that any size dataset will run?
  
  
  
  #### loops over each split in rules_split
  df.ap_L<-list()
  for (i in 1:length(rules_split)){   
    quality(rules_split[[i]])<-cbind(quality(rules_split[[i]]), oddsRatio = interestMeasure(rules_split[[i]], measure="oddsRatio", bin.data))
    #getting out the assays ->phenotype relations of interest
    rules_split[[i]]<-subset(rules_split[[i]], subset = lhs %in% assayItems & lift >1 & rhs %in% phenoItems)
    #convert rules to a data frame object for graphing
    #the output files can be read into cytoscape easily
    df.ap_L[[i]]<-as(rules_split[[i]], "data.frame")
    #currently rules is a factor, setting up to remove redundant rules
    df.ap_L[[i]]$rules<-as.character(df.ap_L[[i]]$rules)
  }
  
  
  
  df.ap<-data.frame()
  
  for (i in 1:length(df.ap_L)){
    df.ap<-rbind(df.ap,df.ap_L[[i]])
  }
  
  colnames(df.ap)<-colnames(df.ap_L[[1]])
  
  temp<-sapply(strsplit(df.ap$rules, split="[{}]"), unlist)
  
  
  
  df.ap$lhs<-temp[2,]          
  df.ap$rhs<-temp[4,]
  
    
  df.ap$lhs<-gsub('=1','',df.ap$lhs)          ####getting rid of the =1 so code will run
  df.ap$rhs<-gsub('=1','',df.ap$rhs)
  
  
####Experimental section: still in beta testing#############  
  
 # print("checkpoint5")
  
#  x<-which(colnames(data_stdlifts) %in% colnames(AssayData))
 # y<-which(colnames(data_stdlifts) %in% colnames(PhenoData))
  
#  assay_stdlifts<-data_stdlifts[x] 
#  pheno_stdlifts<-data_stdlifts[y]
  
 
  
 # df.ap_l<-split(df.ap,cut(seq(nrow(df.ap)),6))   #### the 6 is arbitrary, code will run without splitting data here but takes much longer
  #### like 2+ hours compared to 45 min.  Haven't experimented much yet with other numbers, might speed up
  #### code even more.  Again output is the same regardless of splitting or not
  
  
  #for (k in 1:length(df.ap_l)){
    
    ####gene probabilities###3
  #  lhs_prob<-c()
  #  for (i in 1:nrow(df.ap_l[[k]])){
   #   lhs_prob[i]=sum(assay_stdlifts[,which(colnames(assay_stdlifts)==df.ap_l[[k]]$lhs[i])],na.rm=T)/nrow(assay_stdlifts)
  #  }
    
    ######disease probabilities#####
  #  rhs_prob<-c()
  #  for (i in 1:nrow(df.ap_l[[k]])){
   #   rhs_prob[i]=sum(pheno_stdlifts[,which(colnames(pheno_stdlifts)==df.ap_l[[k]]$rhs[i])],na.rm=T)/nrow(pheno_stdlifts)
  #  }
    
  #  df.ap_l[[k]]<-cbind(df.ap_l[[k]],lhs_prob,rhs_prob )
    
    ##initializing columns
  #  lambda<-vector(length=nrow(df.ap_l[[k]]))
 #   mu<-vector(length=nrow(df.ap_l[[k]]))
 #   standard_lift<-vector(length=nrow(df.ap_l[[k]]))
   # lam2<-vector(length=nrow(df.ap_l[[k]]))
   # standard_lift2<-vector(length=nrow(df.ap_l[[k]]))
  #  df.ap_l[[k]]<-cbind(df.ap_l[[k]],lambda,lam2,mu,standard_lift,standard_lift2)
    
    
    ##standard lift calculations
  #  for (i in 1:nrow(df.ap_l[[k]])){
     # df.ap_l[[k]]$lambda[i] <- (max(df.ap_l[[k]]$lhs_prob[i]+df.ap_l[[k]]$rhs_prob[i]-1,1/nrow(data_stdlifts)))/(df.ap_l[[k]]$lhs_prob[i]*df.ap_l[[k]]$rhs_prob[i])
    #  df.ap_l[[k]]$mu[i] <- 1/(max(df.ap_l[[k]]$lhs_prob[i],df.ap_l[[k]]$rhs_prob[i]))
    #  df.ap_l[[k]]$standard_lift[i] <- (df.ap_l[[k]]$lift[i]-df.ap_l[[k]]$lambda[i])/(df.ap_l[[k]]$mu[i]-df.ap_l[[k]]$lambda[i])
     # df.ap_l[[k]]$lam2[i]<-max((df.ap_l[[k]]$lhs_prob[i]+df.ap_l[[k]]$rhs_prob[i]-1)/(df.ap_l[[k]]$lhs_prob[i]*df.ap_l[[k]]$rhs_prob[i]),(4*sup)/((1+sup)^2),sup/(df.ap_l[[k]]$lhs_prob[i]*df.ap_l[[k]]$rhs_prob[i]),0.001/df.ap_l[[k]]$rhs_prob[i])
    #  df.ap_l[[k]]$standard_lift2[i]<-(df.ap_l[[k]]$lift[i]-df.ap_l[[k]]$lam2[i])/(df.ap_l[[k]]$mu[i]-df.ap_l[[k]]$lam2[i])
   # }
    
  #}
  
  #df.ap2<-data.frame()
  
  #for (i in 1:length(df.ap_l)){
  #  df.ap2<-rbind(df.ap2,df.ap_l[[i]])
 # }
  
  ######################
  #to view your top rules:
  #inspect(head(sort(rules.ap, by = "lift"), n = 10))
  #
  
#####End of experimental section####### 
  
  
  print("------------------------------------------------------------------------------------")
  sink()
  header <- as.data.frame((df.ap))
  header$source <- sapply(strsplit(as.character(header[,8]),'[.]'), "[", 1)
  header$rhs <- sapply(strsplit(as.character(header[,8]),'[.]'), "[", 2)
  
  return(header)
  
}
date <- Sys.Date()
today <- format(date, format="_%B_%d_%Y")

### logfile => this is the location of the logfile where some summary and QC stats will be saved
logfile1<-paste(Root_location,"/OutputFiles/TCgenes_in_CTD_ALLchem_logfile",today,".txt",sep="")
logfile2<-paste(Root_location,"/OutputFiles/allCTD_ALLchem_logfile",today,".txt",sep="")
logfile3<-paste(Root_location,"/OutputFiles/Toxcast_ALLchem_logfile",today,".txt",sep="")

### This is the section that calls the rules and arules function
### The parameters that can be changed in the arules function in line 199
### arules_function <- function(AssayData,PhenoData,sup,minactivity,logfile, datamerge)
### AssayData=> gene data (Toxcast or CTD), 
### PhenoData=> disease data (CTD), 
### sup=> is the minimum support, 
### minactivity=> is the minimum level of activity per row (i.e chemical) or column (i.e. disease or gene)
### logfile => this is the logfile from above

rules1<-arules_function(Gene_data_CTD_TC_human,Disease_data_CTD,10,1,logfile1) #,"CTDTC_ALLchem_merge.txt")
rules2<-arules_function(CTD_all,Disease_data_CTD,10,2,logfile2) #, "CTDALL_ALLchem_merge.txt")
rules3<-arules_function(Gene_data_Toxcast_human,Disease_data_CTD,10,1,logfile3) #,"Toxcast_ALLchem_merge.txt")


#####Writing out to file just the association rules the disease codes (IDs) are used here
write.table(rules1,paste(Root_location,"/OutputFiles/workflow4_TCgenes_in_CTD_ALLchem_assoc_rules",today,".txt",sep=''),sep="\t",row.names=F)
write.table(rules2,paste(Root_location,"/OutputFiles/workflow4_allCTD_ALLchem_assoc_rules",today,".txt",sep=''),sep="\t",row.names=F)
write.table(rules3,paste(Root_location,"/OutputFiles/workflow4_Toxcast_ALLchem_assoc_rules",today,".txt",sep=''),sep="\t",row.names=F)


######Merging disease names with modified result file and Mapping the disease names 
######to the association rules output and Writing out to file
Disease_map <- read.delim(paste(Root_location,"/InputFiles/CTD_disease_map.txt", sep=''),sep='\t', header=TRUE, comment.char='', strip.white=TRUE)
Toxcast_Results_wtDisease_map <- unique(merge(rules3, Disease_map, by.x='rhs', by.y='DiseaseID'))
write.table(Toxcast_Results_wtDisease_map,paste(Root_location,"/OutputFiles/workflow4_Toxcast_ALLchem_assoc_rules_wtDiseaseMap.txt",sep=''), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
CTD_TC_Results_wtDisease_map <- unique(merge(rules1, Disease_map, by.x='rhs', by.y='DiseaseID'))
write.table(CTD_TC_Results_wtDisease_map,paste(Root_location,"/OutputFiles/workflow4_TCgenes_in_CTD_ALLchem_assoc_rules_wtDiseaseMap.txt",sep=''), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
CTD_all_Results_wtDisease_map <- unique(merge(rules2, Disease_map, by.x='rhs', by.y='DiseaseID'))
write.table(CTD_all_Results_wtDisease_map,paste(Root_location,"/OutputFiles/workflow4_allCTD_ALLchem_assoc_rules_wtDiseaseMap.txt",sep=''), sep='\t', col.names=TRUE, row.names=FALSE,quote=FALSE)
