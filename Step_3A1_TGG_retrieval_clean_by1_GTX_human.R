###Author: Noffisat Oki
###Description:
### - Data retrival from TG Gates on a set of pre-selected chemicals.
###   In this particular example we extract human in vitro data for 10 GTX.

#install.packages(magrittr)
#install.packages('jsonlite')
#install.packages(httr)
#install.packages(data.table)
#install.packages(tidyr)
#install.packages(RCurl)
#install.packages(XML)
#install.packages(stringi)
library(magrittr)
library(jsonlite)
library(httr)
set_config(config(ssl_verifypeer = 0L))
library(data.table)
library(tidyr)
library(RCurl)
library(XML)
library(stringi)
## When run separately and locally set your home root location
# Root_location <- "/home/dimiter/Desktop/Advance"

###Getting all the compounds
compounds_disc <- GET("http://open-tggates-api.cloud.douglasconnect.com/v2/compounds?limit=none")
compoundslist <- content(compounds_disc, as= "text", encoding = "UTF-8")
parsed_compounds <- compoundslist%>% fromJSON
compounds <- parsed_compounds$compounds


###Adding Casrn to TGG compounds
caslist <- mat.or.vec(nrow(compounds), 4)#, data=0)
caslist[,c(1:2)]<- as.matrix(compounds[,])
for (row in 1:nrow(caslist)) {
  call2 <- paste('https://cactus.nci.nih.gov/chemical/structure/',caslist[row, 2],'/cas', sep="")
  cas_get <- GET(URLencode(call2))#, reserved = FALSE, repeated = FALSE)
  casrn <- content(cas_get, as= "text", type = "text/html", encoding = "UTF-8")
  parsetest <- htmlParse(casrn, asText = TRUE)
  plain.text <- xpathSApply(parsetest, "//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)][not(ancestor::form)]", xmlValue)#casrn <- parsed_casrn$
  caslist[row,3] <- plain.text #casrn #parsed_casrn <- casrn%>% fromJSON
}

###Adding inchi to TGG
for (row in 1:nrow(caslist)) {
  call3 <- paste('https://cactus.nci.nih.gov/chemical/structure/',caslist[row, 2],'/stdinchi', sep="")
  inchi_get <- GET(URLencode(call3), reserved = FALSE, repeated = FALSE)
  inchi <- content(inchi_get, as= "text", encoding = "UTF-8")
  parsetest2 <- htmlParse(inchi, asText = TRUE)
  plain.text2 <- xpathSApply(parsetest2, "//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)][not(ancestor::form)]", xmlValue)#casrn <- parsed_casrn$
  caslist[row,4] <- plain.text2 #inchi#parsed_casrn <- casrn%>% fromJSON
}
colnames(caslist) <- c("TGG_ID", "TGG_compoundName", "Casrn", "Inchi")
###splitting rows with multiple casrn for same compound
caslist2 <- separate_rows(as.data.frame(caslist), Casrn, sep='\n')
caslist2$Casrn <- as.factor(caslist2$Casrn)
caslist2 <- caslist2[,-4]


###Getting sampleIDs for all samples in TGG
discover_samples <- GET("https://open-tggates-api.cloud.douglasconnect.com/v2/samples?limit=none")
disc_result <- content(discover_samples, as= "text", encoding = "UTF-8")
parsed_disc <- disc_result%>% fromJSON
Samples <- parsed_disc$samples
select_samples <- Samples[,c("sampleId","compoundName", "timepointHr")]
#compoundslist <- unique(select_samples[,2])



###Getting HCC compounds present in TGG. 
### Adjust accordingly to your local directory
HCC_chem_data <- read.table(paste(Root_location,"/InputFiles/Genotoxic_carcinogens_Advance.txt",sep=''),sep='\t',header=T,comment.char='',strip.white=T)
newcat <- which(caslist2$Casrn %in% HCC_chem_data$Cas.number)
newcat2 <- caslist2[newcat,]
TGG_in_HCC_list <- unique(merge(caslist2, HCC_chem_data,  by.x="Casrn", by.y="Cas.number" ,all.x=FALSE, all.y=FALSE))
TGGinHCC_compoundset <- "/OutputFiles/TGG_in_HCC_compounds.txt"
write.table(TGG_in_HCC_list, paste(Root_location,TGGinHCC_compoundset,sep=''),sep="\t",row.names=F)




###############start from here


###Selecting the sample IDs for the HCC compounds in TGG
HCC_samples_in_TGG <- unique(merge(newcat2, Samples, by.x ="TGG_compoundName", by.y = "compoundName" ))
HCC_in_TGG_sampleset <- "/OutputFiles/TGG_in_HCC_samples_GTX.txt"
write.table(HCC_samples_in_TGG[,c(1:11)], paste(Root_location,HCC_in_TGG_sampleset,sep=''),sep="\t",row.names=F)

sampleIDs <- unique(HCC_samples_in_TGG$sampleId)
sampleIDs_human <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Human')),]$sampleId  #&& HCC_samples_in_TGG$tissue =='Liver')),]
sampleIDs_rat <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Rat')),]$sampleId
#sampleIDs.rat.liv.invitro <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Rat' & HCC_samples_in_TGG$tissue =='Liver' & HCC_samples_in_TGG$cellType =='in vitro')),]$sampleId   # && HCC_samples_in_TGG$cellType =='in vivo')),]
#sampleIDs.rat.liv.invivo <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Rat' & HCC_samples_in_TGG$tissue =='Liver' & HCC_samples_in_TGG$cellType =='in vivo')),]$sampleId      # && HCC_samples_in_TGG$cellType =='in vivo')),]
#sampleIDs.rat.kid.invitro <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Rat' & HCC_samples_in_TGG$tissue =='Kidney' & HCC_samples_in_TGG$cellType =='in vitro')),]$sampleId   # && HCC_samples_in_TGG$cellType =='in vivo')),]
#sampleIDs.rat.kid.invivo <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Rat' & HCC_samples_in_TGG$tissue =='Kidney' & HCC_samples_in_TGG$cellType =='in vivo')),]$sampleId     # && HCC_samples_in_TGG$cellType =='in vivo')),]


###Create a loop to paste sample IDs for compounds in HCC set  
#####updates 6/4/2018
###HUMAN###
valuetype1<- "&valueTypeFilter=log2fold|pValue" #"&valueTypeFilter=absolute" 
valuetype2<- "&valueTypeFilter=absolute"                                                                                                     
limit<- "&limit=none" #"&limit=2000000" #

urlhuman1 <- paste("https://open-tggates-api.cloud.douglasconnect.com/v2/results?sampleIdFilter=",sampleIDs_human[1],valuetype1,limit, sep="")
#start_timeA <- timestamp()
output <-  GET(urlhuman1) 
text_result <- content(output, as= "text", encoding = "UTF-8")
#testjson<- fromJSON(text_result, flatten = TRUE)
#testjson <- as.data.frame(testjson)
#plain_results_log2fold_human <- testjson[,c(2:7)]#plain_results#[which(plain_results$valueType =='log2fold'),]
parsed_result <- text_result %>% fromJSON
plain_results <- parsed_result$results
plain_results_log2fold_human <- plain_results#[which(plain_results2$valueType =='absolute'),]

#stop_timeB <- timestamp()

urlhuman2 <- paste("https://open-tggates-api.cloud.douglasconnect.com/v2/results?sampleIdFilter=",sampleIDs_human[1],valuetype2,limit, sep="")
output2 <-  GET(urlhuman2) 
text_result2 <- content(output2, as= "text", encoding = "UTF-8")
parsed_result2 <- text_result2 %>% fromJSON
plain_results2 <- parsed_result2$results
plain_results_absolute_human <- plain_results2#[which(plain_results2$valueType =='absolute'),]

start_time <- timestamp()
counts=2
print('human counts')
while (counts <= length(sampleIDs_human)){#IDs_initial) {
  urlhuman <- paste("https://open-tggates-api.cloud.douglasconnect.com/v2/results?sampleIdFilter=",sampleIDs_human[counts], sep="")
  #print("counts")
  print(counts)
  
  #start_time1 <- timestamp()
  urlhuman1 <- paste(urlhuman,valuetype1,limit, sep="")
  output <-  GET(urlhuman1) 
  text_result <- content(output, as= "text", encoding = "UTF-8")
  #testjson2<- fromJSON(text_result, flatten = TRUE)
  #testjson2 <- as.data.frame(testjson2)
  #plain_results_log2fold_human <- rbind(plain_results_log2fold_human, testjson2[,c(2:7)])
  #stop_time1 <- timestamp()  
  parsed_result <- text_result%>% fromJSON
  b <- parsed_result$results
  plain_results_log2fold_human <- rbind(plain_results_log2fold_human, b)
  
  #start_time2 <- timestamp()
  urlhuman2 <- paste(urlhuman,valuetype2,limit, sep="")
  output2 <-  GET(urlhuman2) 
  text_result2 <- content(output2, as= "text", encoding = "UTF-8")
  parsed_result2 <- text_result2%>% fromJSON
  b <- parsed_result2$results
  plain_results_absolute_human <- rbind(plain_results_absolute_human, b)
  #stop_time2 <- timestamp()
  
  counts <- counts+1 #10
  #print('counts')
  #print(counts)

}

human_dataset1 <- "/OutputFiles/TGG_HCC_data_absolute_human_GTX.txt"
human_dataset2 <- "/OutputFiles/TGG_HCC_data_log2fold_human_GTX.txt"
plain_results_absolute_human$geneSymbols <- unlist(as.character(plain_results_absolute_human$geneSymbols))
write.table(plain_results_absolute_human, paste(Root_location,human_dataset1,sep=''),sep="\t",row.names=F)
plain_results_log2fold_human$geneSymbols <- unlist(as.character(plain_results_log2fold_human$geneSymbols))
write.table(plain_results_log2fold_human, paste(Root_location,human_dataset2,sep=''),sep="\t",row.names=F)

#testjson3 <- testjson2
#testjson3$results.geneSymbols <- unlist(as.character(testjson3$results.geneSymbols))

stoptime <- timestamp()

