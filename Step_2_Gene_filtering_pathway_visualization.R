#Author: Tatyana Doktorova and Noffisat Oki
#Version: 12_12_2018
# The script filters the most important genes in the case of HCC and facilitates their visulaisation
#
####            
#####################################################################

library(shiny)
#install.packages("C:\\Users\\tatyana\\Documents\\Advance project\\shiny_1.1.0.tar.gz", repos = NULL, type="source")
library(igraph)
library (dplyr)
library(VennDiagram)
library(treemap)
library(data.tree)
library(data.table)

library(DiagrammeR)
library(collapsibleTree)
library(RColorBrewer)

setwd("C:\\Users\\tatyana\\Documents\\Projects\\Advance project")

#AWF: If the files are to be run in a workflow this step is not needed: Upload disease relevant files

#TC_perc <- read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/TC_perc.txt",sep = '\t', header=TRUE)
#CTD_perc <- read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/CTD_perc.txt",sep = '\t', header=TRUE)
#CTDTC_perc <- read.csv("file:///C:/Users/tatyana/Documents/Projects/Advance project/CTDTC_perc.txt",sep = '\t', header=TRUE)

#Select the unique categories of interest

TC_perc_sub<- subset(TC_perc, select=c("Parent.Name","lhs" , "Major", "Minor"))
CTD_perc_sub<- subset(CTD_perc, select=c("Parent.Name","lhs" , "Major", "Minor"))
CTDTC_perc_sub<- subset(CTDTC_perc, select=c("Parent.Name","lhs",  "Major", "Minor"))

TC_perc_parent = unique(TC_perc_sub)
CTD_perc_parent = unique(CTD_perc_sub)
CTDTC_perc_parent = unique(CTDTC_perc_sub)

# Indentify unique genes only 

TC_gene<- unique(subset(TC_perc, select=c("lhs"  )))
CTDTC_gene<- unique(subset(CTDTC_perc, select=c("lhs"  )))
CTD_gene<- unique(subset(CTD_perc, select=c("lhs"  )))

TC_perc_gene = TC_perc_sub %>% distinct(lhs)
CTD_perc_gene = CTD_perc_sub %>% distinct(lhs)
CTDTC_perc_gene = CTDTC_perc_sub %>% distinct(lhs)

#Look for overlapping genes

TC_CTDTC_CTD_overlap<- Reduce(intersect, list(TC_perc_gene$lhs,CTD_perc_gene$lhs, CTDTC_perc_gene$lhs ))
TC_CTDTC_overlap<-Reduce(intersect, list(TC_perc_gene$lhs, CTDTC_perc_gene$lhs ))
TC_CTD_overlap<-Reduce(intersect, list(TC_perc_gene$lhs, CTD_perc_gene$lhs ))
CTDTC_CTD_overlap<-Reduce(intersect, list(CTDTC_perc_gene$lhs, CTD_perc_gene$lhs ))

#Merge all unique genes
all<-unique(rbind(TC_gene,CTDTC_gene, CTD_gene))

#Overlaps as dataframes
TC_CTDTC_CTD_overlap_df <- as.data.frame(TC_CTDTC_CTD_overlap)
TC_CTDTC_overlap_df <- as.data.frame(TC_CTDTC_overlap)
TC_CTD_overlap_df <- as.data.frame(TC_CTD_overlap)
CTDTC_CTD_overlap_df <- as.data.frame(CTDTC_CTD_overlap)


#Merge the genes plus respective pathway info together
merge_all<-unique(rbindlist(list(CTD_perc_parent,TC_perc_parent,CTDTC_perc_parent)))
#write.csv(merge_all, file="HCC_Toxcast_CTD.csv")
HCC_Toxcast_CTD<-as.data.frame (merge_all)

ALL = merge_all %>% distinct(lhs)
gene_major_subset<-unique(data.frame(merge_all$"lhs",merge_all$"Major"))


## Visualization: 1 
###Tree diagram visualization starting with Major/Minor/Parent_Name/gene

merge_all$pathString <- paste("All pathways_gene relationships", 
                            merge_all$Major, 
                            merge_all$Minor,   merge_all$Parent.Name,  merge_all$lhs,
                            sep = "/")

tree_diagram <- as.Node(merge_all)
print(tree_diagram)
#export_graph(ToDiagrammeRGraph(tree_diagram), "All_pathways_gene_relationships.pdf")
#write.csv(tree_diagram, file="All_pathways_gene_relationships.csv")
#write.table(tree_diagram, file="All_pathways_gene_relationships.txt", sep = "\t")

### Tree diagram visualization starting with gene/Major category

merge_all$pathString <- paste("Key gene_major pathway relationships", 
                              merge_all$lhs,merge_all$Major, 
                              sep = "/")

tree_diagram <- as.Node(merge_all)
print(tree_diagram)

#write.csv(tree_diagram, file="Key_gene_major_pathway relationships.csv")
#write.table(tree_diagram, file="Key_gene_major_pathway relationships.txt", sep = "\t")
#export_graph(ToDiagrammeRGraph(tree_diagram), "Key_gene_major_pathway_relationships.pdf")

 ##Visualization 2: Collapsible tree (interactive)
### Visualization of all HCC data
HCC<-data.frame(merge_all)
collapsibleTree(df=HCC, c("Major", "Minor","Parent.Name", "lhs"),  fill = "lightgreen")
collapsibleTree(df=merge_all, c("Major",  "lhs"),  fill = "lightpink")             
collapsibleTree(df=merge_all, c("lhs","Parent.Name","Minor","Major"),  fill = "green")

###Visualization of individual genes within the HCC subset
GENE<- filter(merge_all,lhs == "ACTB")
collapsibleTree(df=GENE, c("lhs","Major","Minor","Parent.Name"),  fill = "green")

GENE1<- filter(merge_all,lhs == "XRCC1")
collapsibleTree(df=GENE1, c("lhs","Major","Minor","Parent.Name"),  fill = "green")
