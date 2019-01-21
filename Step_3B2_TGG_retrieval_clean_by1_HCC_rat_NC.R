###Author: Noffisat Oki
###Description:
### - Data retrival from TG Gates on a set of pre-selected chemicals.
###   In this particular example we extract human in vitro data for 10 NC.



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
Root_location <- "/home/dimiter/Desktop/Advance"

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



###Getting HCC compounds present in TGG
HCC_chem_data <- read.table(paste(Root_location,"/InputFiles/Non_carcinogens_Advance.txt",sep=''),sep='\t',header=T,comment.char='',strip.white=T)
newcat <- which(caslist2$Casrn %in% HCC_chem_data$Cas.number)
newcat2 <- caslist2[newcat,]
TGG_in_HCC_list <- unique(merge(caslist2, HCC_chem_data,  by.x="Casrn", by.y="Cas.number" ,all.x=FALSE, all.y=FALSE))
TGGinHCC_compoundset <- "/OutputFiles/TGG_in_HCC_compounds.txt"
write.table(TGG_in_HCC_list, paste(Root_location,TGGinHCC_compoundset,sep=''),sep="\t",row.names=F)



###############start from here


###Selecting the sample IDs for the HCC compounds in TGG
HCC_samples_in_TGG <- unique(merge(newcat2, Samples, by.x ="TGG_compoundName", by.y = "compoundName" ))
HCC_in_TGG_sampleset <- "/OutputFiles/TGG_in_HCC_samples_rat_NC.txt"
write.table(HCC_samples_in_TGG[,c(1:11)], paste(Root_location,HCC_in_TGG_sampleset,sep=''),sep="\t",row.names=F)

sampleIDs <- unique(HCC_samples_in_TGG$sampleId)
sampleIDs_human <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Human')),]$sampleId  #&& HCC_samples_in_TGG$tissue =='Liver')),]
sampleIDs_rat <- HCC_samples_in_TGG[c(which(HCC_samples_in_TGG$organism =='Rat')),]$sampleId


###RAT###
valuetype1<- "&valueTypeFilter=log2fold|pValue" #"&valueTypeFilter=absolute" 
#valuetype2<- "&valueTypeFilter=absolute" 
limit<- "&limit=none" 
urlrat1 <- paste("https://open-tggates-api.cloud.douglasconnect.com/v2/results?sampleIdFilter=",sampleIDs_rat[1],valuetype1,limit, sep="")
output <-  GET(urlrat1) 
text_result <- content(output, as= "text", encoding = "UTF-8")
parsed_result <- text_result %>% fromJSON
plain_results <- parsed_result$results
plain_results_log2fold_rat <- plain_results#[which(plain_results$valueType =='log2fold'),]

#urlrat2 <- paste("https://open-tggates-api.cloud.douglasconnect.com/v2/results?sampleIdFilter=",sampleIDs_rat[1],valuetype2,limit, sep="")
#output2 <-  GET(urlrat2) 
#text_result2 <- content(output2, as= "text", encoding = "UTF-8")
#parsed_result2 <- text_result2 %>% fromJSON
#plain_results2 <- parsed_result2$results
#plain_results_absolute_rat <- plain_results2#[which(plain_results2$valueType =='absolute'),]

start_time <- timestamp()
counts=2
print('rat counts')
#print(counts)
#IDs_leftover <- length(sampleIDs_rat)%%10
#IDs_initial <- length(sampleIDs_rat)- IDs_leftover
while (counts <= length(sampleIDs_rat)){#IDs_initial) {
  urlrat <- paste("https://open-tggates-api.cloud.douglasconnect.com/v2/results?sampleIdFilter=",sampleIDs_rat[counts], sep="")
  #print("counts")
  print(counts)
  
  urlrat1 <- paste(urlrat,valuetype1,limit, sep="")
  output <-  GET(urlrat1) 
  text_result <- content(output, as= "text", encoding = "UTF-8")
  parsed_result <- text_result%>% fromJSON
  b <- parsed_result$results
  plain_results_log2fold_rat <- rbind(plain_results_log2fold_rat, b)
  
  #urlrat2 <- paste(urlrat,valuetype2,limit, sep="")
  #output2 <-  GET(urlrat2) 
  #text_result2 <- content(output2, as= "text", encoding = "UTF-8")
  #parsed_result2 <- text_result2%>% fromJSON
  #b <- parsed_result2$results
  # plain_results_absolute_rat <- rbind(plain_results_absolute_rat, b)
  counts <- counts+1
  #print('counts')
  #print(counts)
}

#rat_dataset1 <- "/OutputFiles/TGG_data_absolute_rat_HCC.txt"
rat_dataset2 <- "/OutputFiles/TGG_data_log2fold_HCC_rat_NC.txt"
#plain_results_absolute_rat$geneSymbols <- unlist(as.character(plain_results_absolute_rat$geneSymbols))
plain_results_log2fold_rat$geneSymbols <- unlist(as.character(plain_results_log2fold_rat$geneSymbols))
#write.table(plain_results_absolute_rat, paste(Root_location,rat_dataset1,sep=''),sep="\t",row.names=F)
write.table(plain_results_log2fold_rat, paste(Root_location,rat_dataset2,sep=''),sep="\t",row.names=F)
stop_time <- timestamp()


###########################################################

