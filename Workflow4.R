####################################################################
####Author: Noffisat Oki
####Version: 1.0 9-10-2015
####Description: Association rule analysis workflow.
####Notes: Current version analyzes only human data from both the CTD and ToxCast datasets. 
####       Input data and preceding scripts may be modified for analysis of other desired species 
####       present in the ToxCast data.
####       The path to the 'cpAOP_paper' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 17
####       should be used for this purpose.
####Potential issues:
#####################################################################

library(arules)

### Root folder location 
Root_location <- "/share/home/noki/cpAOP_paper_submission"  ##The root folder should be specified here

#Output directory
outDir <- paste(Root_location, "/cpAOP_paper/OutputFiles/", sep='')

######Toxcast
Gene_data_Toxcast_human <- read.csv(paste(Root_location,"/cpAOP_paper/IntermediateData/ToxCast_chem_x_gene_Hitcalls_wt_all_IDs_human.csv",sep=''), sep=',', header=TRUE, comment.char='', strip.white=TRUE)
Gene_data_Toxcast_human <- Gene_data_Toxcast_human[,c(-1,-3,-4,-5)]
###Remove chemicals with no CID (CID=NA) as the gene and disease datasets are merged on CID
Gene_data_Toxcast_human <- na.omit(Gene_data_Toxcast_human)
ids <- Gene_data_Toxcast_human[,1]
##making sure all values are discretized either 1 or 0
Gene_data_Toxcast_human[Gene_data_Toxcast_human>=1] <- 1
Gene_data_Toxcast_human[Gene_data_Toxcast_human<1] <- 0
Gene_data_Toxcast_human <- Gene_data_Toxcast_human[,c(-1)]
Gene_data_Toxcast_human <- cbind(ids, Gene_data_Toxcast_human)
##Splitting the data to run in smaller batches due to memory limitations
Gene_data_Toxcast_human1 <- Gene_data_Toxcast_human[,1:(ncol(Gene_data_Toxcast_human)/2)]
Gene_data_Toxcast_human2 <- Gene_data_Toxcast_human[,((ncol(Gene_data_Toxcast_human)/2)+1):ncol(Gene_data_Toxcast_human)]
Gene_data_Toxcast_human2 <- cbind(Gene_data_Toxcast_human$ids,Gene_data_Toxcast_human2 )
colnames(Gene_data_Toxcast_human1)[1] <- "cid"
colnames(Gene_data_Toxcast_human2)[1] <- "cid"

######CTD#####
###CTD genes in Toxcast
Gene_data_CTD_TC_human <- read.csv(paste(Root_location, "/cpAOP_paper/IntermediateData/CTD_chem_x_gene_human_wt_all_IDs_humanTCgenes.csv",sep=''), sep=',', header=TRUE, comment.char='', strip.white=TRUE)
Gene_data_CTD_TC_human <- na.omit(Gene_data_CTD_TC_human)
cid <- Gene_data_CTD_TC_human[,3]
Gene_data_CTD_TC_human <- Gene_data_CTD_TC_human[,c(-1,-2,-3)]
##making sure all values are discretized either 1 or 0
Gene_data_CTD_TC_human[Gene_data_CTD_TC_human>=1] <- 1
Gene_data_CTD_TC_human[Gene_data_CTD_TC_human<1] <- 0
Gene_data_CTD_TC_human <- cbind(cid, Gene_data_CTD_TC_human)
Gene_data_CTD_TC_human1 <- Gene_data_CTD_TC_human[,1:(ncol(Gene_data_CTD_TC_human)/2)]
Gene_data_CTD_TC_human2 <- Gene_data_CTD_TC_human[,((ncol(Gene_data_CTD_TC_human)/2)+1):ncol(Gene_data_CTD_TC_human)]
Gene_data_CTD_TC_human2 <- cbind(Gene_data_CTD_TC_human$cid,Gene_data_CTD_TC_human2 )
colnames(Gene_data_CTD_TC_human2)[1] <- "cid"

###All diseases in CTD
Disease_data_CTD <- read.csv(paste(Root_location, "/cpAOP_paper/IntermediateData/CTD_chem_x_disease_wt_all_IDs.csv",sep=''), sep=',', header=TRUE, comment.char='', strip.white=TRUE)
Disease_data_CTD <- na.omit(Disease_data_CTD)
Disease_data_CTD <- Disease_data_CTD[,c(-1,-2)]
cids <- Disease_data_CTD[,1]
Disease_data_CTD <- Disease_data_CTD[,c(-1)]
##making sure all values are discretized either 1 or 0
Disease_data_CTD[Disease_data_CTD>=1] <- 1
Disease_data_CTD[Disease_data_CTD<1] <- 0
Disease_data_CTD <- cbind(cids, Disease_data_CTD)
##Splitting the data to run in smaller batches due to memory limitations
Disease_data_CTD_split1 <- Disease_data_CTD[,1:(ncol(Disease_data_CTD)/4)]
Disease_data_CTD_split2 <- Disease_data_CTD[,((ncol(Disease_data_CTD)/4)+1):(ncol(Disease_data_CTD)/2)]
Disease_data_CTD_split3 <- Disease_data_CTD[,((ncol(Disease_data_CTD)/2)+1):((ncol(Disease_data_CTD)/2)+(ncol(Disease_data_CTD)/4))]
Disease_data_CTD_split4 <- Disease_data_CTD[,(((ncol(Disease_data_CTD)/2)+(ncol(Disease_data_CTD)/4))+1):ncol(Disease_data_CTD)]
Disease_data_CTD_split2 <- cbind(Disease_data_CTD$cid,Disease_data_CTD_split2 )
Disease_data_CTD_split3 <- cbind(Disease_data_CTD$cid,Disease_data_CTD_split3 )
Disease_data_CTD_split4 <- cbind(Disease_data_CTD$cid,Disease_data_CTD_split4 )
colnames(Disease_data_CTD_split1)[1] <- "cid"
colnames(Disease_data_CTD_split2)[1] <- "cid"
colnames(Disease_data_CTD_split3)[1] <- "cid"
colnames(Disease_data_CTD_split4)[1] <- "cid"

####Output file names
Toxcast_gene_diseases_phenorules <- "Toxcast_genedisease_human_PhenoRulesGraph"
CTD_gene_disease_TCgenes_phenorules <- "CTD_genedisease_human_PhenoRulesGraph_TCgenes"

###association mining function
arules_function <- function(AssayData,PhenoData,sup,logfile){
  
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
  dataBinS<-dataBin[,!apply(dataBin,2,function(x) sum(abs(x), na.rm=TRUE) <= 2)]
  #remove rows which are all '0' indicating no effect for a chemical
  dataBinS<-dataBin[!apply(dataBin,1,function(x) sum(abs(x), na.rm=TRUE)<= 2),]
   

  ###Adding these two extra lines to eliminate the possibility of inf when calculating odds ratio
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
  print("minimum support %)")
  sup <- sup/nrow(dataBinS)
  

  #to simplify the rules, removing the level '0' and replacing it with NA
  #the effect will be to have just presence for a 'response' measure
  dataBinEffect<-dataBinS; dataBinS<-as.data.frame(apply(dataBinS,2,as.factor))
  dataBinEffect[dataBinEffect==0]<-NA; dataBinEffect<-as.data.frame(apply(dataBinEffect,2,as.factor))
  #create the transaction object
  bin.data<-as(dataBinEffect, "transactions")
  #this gets out the "labels" for the rules coming from each dataset
  assayItems<-intersect(colnames(bin.data), paste(assays,1, sep='='))
  phenoItems<-intersect(colnames(bin.data), paste(pheno,1, sep='='))
  
  #############################################
  #Support = #instance of your 'set'/ total # transactions
  #Confidence = support(x,y)/support(LHS)
  #Lift = support(x,y)/support(x)*support(y) (ratio of observed support to that if X and Y were independent)
  #setting the min/max len=2 will limit to 2 itemset rules
  rules<-apriori(bin.data, parameter =list(support = sup, confidence = 0.001, minlen=2, maxlen=2))
  rules.ap<-subset(rules, subset = lift >1)
  quality(rules.ap)<-cbind(quality(rules.ap), oddsRatio = interestMeasure(rules.ap, method="oddsRatio", bin.data))
  #getting out the assays ->phenotype relations of interest
  rules.ap<-subset(rules.ap, subset = lhs %in% assayItems & lift >1 & rhs %in% phenoItems)
  #convert rules to a data frame object for graphing
  #the output files can be read into cytoscape easily
  df.ap<-as(rules.ap, "data.frame")
  #currently rules is a factor, setting up to remove redundant rules
  df.ap$rules<-as.character(df.ap$rules)
  temp<-sapply(strsplit(df.ap$rules, split="[{}]"), unlist)
  df.ap$lhs<-temp[2,]
  df.ap$rhs<-temp[4,]

  ######################
  #to view your top rules:
  #inspect(head(sort(rules.ap, by = "lift"), n = 10))
  #


  print("------------------------------------------------------------------------------------")
  sink()
  return(df.ap)
}



################ToxCast data and CTD disease 4 splits of disease data 
date <- Sys.Date()
today <- format(date, format="_%B_%d_%Y")
logfile1 <- paste(outDir,"Gene_disease_Toxcast_human_all_logfile_split4",today, ".txt", sep="")
Gene_disease_Toxcast_human1a_4split <- arules_function(Gene_data_Toxcast_human1, Disease_data_CTD_split1,10,logfile1)
Gene_disease_Toxcast_human1b_4split <- arules_function(Gene_data_Toxcast_human1, Disease_data_CTD_split2,10,logfile1)
Gene_disease_Toxcast_human1c_4split <- arules_function(Gene_data_Toxcast_human1, Disease_data_CTD_split3,10,logfile1)
Gene_disease_Toxcast_human1d_4split <- arules_function(Gene_data_Toxcast_human1, Disease_data_CTD_split4,10,logfile1)
Gene_disease_Toxcast_human2a_4split <- arules_function(Gene_data_Toxcast_human2, Disease_data_CTD_split1,10,logfile1)
Gene_disease_Toxcast_human2b_4split <- arules_function(Gene_data_Toxcast_human2, Disease_data_CTD_split2,10,logfile1)
Gene_disease_Toxcast_human2c_4split <- arules_function(Gene_data_Toxcast_human2, Disease_data_CTD_split3,10,logfile1)
Gene_disease_Toxcast_human2d_4split <- arules_function(Gene_data_Toxcast_human2, Disease_data_CTD_split4,10,logfile1)
##combining all function return outputs
Gene_disease_Toxcast_human_all_split <- rbind(Gene_disease_Toxcast_human1a_4split,Gene_disease_Toxcast_human1b_4split,Gene_disease_Toxcast_human1c_4split,Gene_disease_Toxcast_human1d_4split,
                                        Gene_disease_Toxcast_human2a_4split,Gene_disease_Toxcast_human2b_4split, Gene_disease_Toxcast_human2c_4split,Gene_disease_Toxcast_human2d_4split)
write.table(Gene_disease_Toxcast_human_all_split,file=file.path(outDir, paste(Toxcast_gene_diseases_phenorules,today,"_split4.csv",sep="")), sep=',', col.names=TRUE,row.names=FALSE,quote=FALSE )

###########CTD data ToxCast gene subset and CTD disease 4 splits of disease data 
date <- Sys.Date()
today <- format(date, format="_%B_%d_%Y")
logfile2 <- paste(outDir,"Gene_disease_CTD_TCgenes_human_all_logfile_split4",today,".txt", sep="")
Gene_disease_CTD_TCgenes_human1a_4split <- arules_function(Gene_data_CTD_TC_human1, Disease_data_CTD_split1,10,logfile2)
Gene_disease_CTD_TCgenes_human1b_4split <- arules_function(Gene_data_CTD_TC_human1, Disease_data_CTD_split2,10,logfile2)
Gene_disease_CTD_TCgenes_human1c_4split <- arules_function(Gene_data_CTD_TC_human1, Disease_data_CTD_split3,10,logfile2)
Gene_disease_CTD_TCgenes_human1d_4split <- arules_function(Gene_data_CTD_TC_human1, Disease_data_CTD_split4,10,logfile2)
Gene_disease_CTD_TCgenes_human2a_4split <- arules_function(Gene_data_CTD_TC_human2, Disease_data_CTD_split1,10,logfile2)
Gene_disease_CTD_TCgenes_human2b_4split <- arules_function(Gene_data_CTD_TC_human2, Disease_data_CTD_split2,10,logfile2)
Gene_disease_CTD_TCgenes_human2c_4split <- arules_function(Gene_data_CTD_TC_human2, Disease_data_CTD_split3,10,logfile2)
Gene_disease_CTD_TCgenes_human2d_4split <- arules_function(Gene_data_CTD_TC_human2, Disease_data_CTD_split4,10,logfile2)
##combining all function return outputs
Gene_disease_CTD_TCgenes_human_all_split <- rbind(Gene_disease_CTD_TCgenes_human1a_4split,Gene_disease_CTD_TCgenes_human1b_4split,Gene_disease_CTD_TCgenes_human1c_4split,Gene_disease_CTD_TCgenes_human1d_4split, 
                                            Gene_disease_CTD_TCgenes_human2a_4split,Gene_disease_CTD_TCgenes_human2b_4split,Gene_disease_CTD_TCgenes_human2c_4split,Gene_disease_CTD_TCgenes_human2d_4split)
CTD_gene_disease_TCgenes_phenorules_all <- "CTD_genedisease_human_PhenoRulesGraph_TCgenes_all"
write.table(Gene_disease_CTD_TCgenes_human_all_split,file=file.path(outDir, paste(CTD_gene_disease_TCgenes_phenorules_all,today,"_split4.csv",sep="")), sep=',', col.names=TRUE,row.names=FALSE,quote=FALSE )


